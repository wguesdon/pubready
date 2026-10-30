"""Deterministic, non-statistical advice from a versioned decision tree."""

import hashlib
import json
import os
from pathlib import Path
import shlex

import numpy as np
import pandas as pd
import yaml

from .io import read_tidy
from .run import inspect_column, inspect_id_candidate, inspect_labels
from .spec import PLOT_VALUE_OPTIONS, check_recipe_options, read_recipe_options, split_list_option
from .util import level_order, time_level_order

REFERENCE = Path(os.environ.get("PUBREADY_REFERENCE", "/opt/pubready/reference"))
STATUS_CODES = {"complete": 0, "usage_error": 2, "needs_answers": 3,
                "conflict": 4, "out_of_scope": 5}
COLUMN_FLAGS = ("x", "y", "fill", "facet", "id", "group", "time", "event", "label",
                "ci_low", "ci_high", "se", "covariates", "features", "markers")
LIST_COLUMNS = {"covariates", "features", "markers"}
OPERATORS = {"equals", "at_least", "at_most", "subset_of", "not_subset_of", "empty"}


class AdviceError(ValueError):
    """An invalid command or contradictory answer.

    Attributes:
        status: The public result status.
        detail: Structured information about the error.
    """

    def __init__(self, message, status="usage_error", **detail):
        """Initialize an advice error.

        Args:
            message: A complete explanation.
            status: The result status.
            **detail: Context for the conflict record.
        """
        super().__init__(message)
        self.status = status
        self.detail = detail


def load_tree(path=None):
    """Read and validate the tree.

    Args:
        path: Optional YAML path.

    Returns:
        tuple: The tree and its SHA256 digest.
    """
    raw = Path(path or REFERENCE / "decision_tree.yaml").read_bytes()
    tree = yaml.safe_load(raw)
    validate_tree(tree)
    return tree, hashlib.sha256(raw).hexdigest()


def validate_tree(tree):
    """Check the schema, branch coverage, references, reachability and cycles.

    Args:
        tree: The parsed tree.

    Raises:
        ValueError: If the tree is not supported or internally consistent.
    """
    if not isinstance(tree, dict) or tree.get("schema_version") != 1:
        raise ValueError("The advise tree has an unsupported schema version.")
    questions, nodes = tree["questions"], tree["nodes"]
    for key, q in questions.items():
        if not isinstance(key, str) or q["source"] not in ("data", "scientist", "data_or_scientist"):
            raise ValueError("The advise tree has an invalid question.")
        if any(not isinstance(a, str) for a in q.get("answers", {})):
            raise ValueError("Every answer key must be a string.")
        for rule in q.get("data_rules", []):
            if len(OPERATORS.intersection(rule)) != 1:
                raise ValueError("A data rule needs exactly one supported operator.")
    table = read_recipe_options()
    for key, node in nodes.items():
        if not isinstance(key, str):
            raise ValueError("Every node key must be a string.")
        if "ask" in node:
            if node["ask"] not in questions:
                raise ValueError(f"Unknown question in node {key}.")
            if set(node["branches"]) != set(questions[node["ask"]]["answers"]):
                raise ValueError(f"Incomplete answer coverage in node {key}.")
        elif "leaf" in node:
            leaf = node["leaf"]
            if leaf["recipe"] not in table:
                raise ValueError(f"Unknown recipe in node {key}.")
            check_recipe_options(leaf["recipe"], leaf["flags"]["required"])
            check_recipe_options(leaf["recipe"], leaf.get("fixed_flags", {}))
            if any(q not in questions for q in leaf.get("ask_at_leaf", [])):
                raise ValueError(f"Unknown leaf question in {key}.")
        elif "out_of_scope" not in node:
            raise ValueError(f"Invalid node {key}.")
    visited = set()
    pending = [("start", ())]
    while pending:
        key, ancestors = pending.pop()
        if key not in nodes or key in ancestors:
            raise ValueError(f"Missing node or cycle at {key}.")
        visited.add(key)
        for target in nodes[key].get("branches", {}).values():
            pending.append((target, ancestors + (key,)))
    if visited != set(nodes):
        raise ValueError("The advise tree contains unreachable nodes.")


def match_rule(rule, facts):
    """Evaluate one closed-set rule without executing expressions.

    Args:
        rule: A tree rule.
        facts: Facts keyed by dotted names.

    Returns:
        bool: Whether the rule holds.
    """
    value = facts.get(rule["fact"])
    if value is None:
        return False
    if "equals" in rule:
        return value == rule["equals"]
    if "at_least" in rule:
        return value >= rule["at_least"]
    if "at_most" in rule:
        return value <= rule["at_most"]
    if "empty" in rule:
        return (not bool(value)) == rule["empty"]
    if "subset_of" in rule:
        return bool(value) and set(value).issubset(rule["subset_of"])
    return bool(value) and not set(value).issubset(rule["not_subset_of"])


def add_answer(answers, key, value):
    """Add an answer, rejecting contradictory duplicates.

    Args:
        answers: Mutable answer mapping.
        key: Question or plot flag name.
        value: Scalar answer.
    """
    if isinstance(value, (dict, list)) or value is None:
        raise AdviceError(f"Answer {key} must be a scalar value.")
    value = str(value).lower() if isinstance(value, bool) else str(value)
    if key in answers and answers[key] != value:
        raise AdviceError(f"Question {key} has two different answers.", question=key)
    answers[key] = value


class AnswerLoader(yaml.SafeLoader):
    """A safe YAML loader that preserves duplicate keys for validation."""


def answer_mapping(loader, node):
    """Construct an answer mapping without silently replacing duplicate keys.

    Args:
        loader: Safe YAML loader.
        node: Mapping node.

    Returns:
        dict: Validated answers.
    """
    result = {}
    for key_node, value_node in node.value:
        add_answer(result, loader.construct_object(key_node), loader.construct_object(value_node))
    return result


AnswerLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, answer_mapping)


def merge_answers(tree, options, entries=(), filename=None):
    """Combine a file, repeated answers and explicit plot flags.

    Args:
        tree: The decision tree.
        options: Explicit plot flags only.
        entries: KEY=VALUE arguments.
        filename: Optional YAML or JSON file.

    Returns:
        dict: The combined answers.
    """
    answers = {}
    if filename:
        with open(filename, encoding="utf-8") as handle:
            loaded = yaml.load(handle, Loader=AnswerLoader)
        if not isinstance(loaded, dict):
            raise AdviceError("The answers file must contain a mapping.")
        for key, value in loaded.items():
            add_answer(answers, key, value)
    for entry in entries:
        if "=" not in entry:
            raise AdviceError("Each --answer needs KEY=VALUE.")
        key, value = entry.split("=", 1)
        add_answer(answers, key, value)
    accepted = set(tree["questions"]) | set().union(*map(set, read_recipe_options().values()))
    for key in answers:
        if key not in accepted or key == "recipe":
            raise AdviceError(f"Unknown answer key {key}.")
    for key, value in options.items():
        if key in tree["questions"] or key in answers:
            add_answer(answers, key, value)
    if options.get("paired") or answers.get("paired") == "true":
        add_answer(answers, "design", "paired")
    for key, value in answers.items():
        q = tree["questions"].get(key, {})
        if "answers" in q and value not in q["answers"]:
            raise AdviceError(f"Invalid {key}={value}. Choose from {', '.join(q['answers'])}.")
    return answers


def column_facts(df, column, role):
    """Extend shared inspect facts for a selected column.

    Args:
        df: Input table.
        column: Column name.
        role: Plot flag name.

    Returns:
        dict: Dotted facts for the role.
    """
    series = df[column]
    info = inspect_column(column, series)
    labels = inspect_labels(series).dropna()
    levels = level_order(labels)
    numeric = info["type"] == "numeric" and bool(np.isfinite(series.dropna().astype(float)).all())
    result = {"column": column, "kind": info["type"], "numeric": numeric,
              "n_unique": info["n_unique"], "n_missing": info["n_missing"],
              "levels": levels, "n_levels": len(levels),
              "n_per_level": {v: int((labels == v).sum()) for v in levels},
              "min_per_level": int(labels.value_counts().min()) if len(labels) else 0,
              "continuous": numeric and info["n_unique"] > 6}
    if numeric:
        result["all_positive"] = bool(len(series.dropna()) and (series.dropna() > 0).all())
    if role == "event":
        codes = ["1" if v == "TRUE" else "0" if v == "FALSE" else v for v in levels]
        result["codes"] = level_order(codes)
        result["nonzero_codes"] = [v for v in result["codes"] if v != "0"]
    return {f"{role}.{key}": value for key, value in result.items()}


def collect_facts(df, options):
    """Reuse inspect functions and add the facts needed for the walk.

    Args:
        df: Input table.
        options: Selected plot flags.

    Returns:
        dict: Facts keyed by dotted names.
    """
    facts = {"columns": list(df.columns), "id.candidates": []}
    for role in COLUMN_FLAGS:
        if not options.get(role):
            continue
        value = options[role]
        columns = split_list_option(value) if role in LIST_COLUMNS or (role == "x" and value not in df) else [value]
        for column in columns:
            if column not in df:
                raise AdviceError(f"Column {column} from --{role} does not exist.")
        if len(columns) == 1:
            facts.update(column_facts(df, columns[0], role))
    x, y = options.get("x"), options.get("y")
    if x in df:
        for col in df:
            if col in (x, y):
                continue
            candidate = inspect_id_candidate(inspect_labels(df[col]), inspect_labels(df[x]), col, x)
            if candidate:
                facts["id.candidates"].append(candidate)
    if facts["id.candidates"]:
        facts["id.proposal"] = "paired" if facts.get("x.n_levels") == 2 else "repeated"
    facts["n_factors"] = sum(bool(options.get(k)) for k in ("x", "fill", "facet"))
    if y and not options.get("facet"):
        facts["layout"] = "one"
    elif not y:
        numeric = [c for c in df if c not in (x, options.get("id"), options.get("label"))
                   and inspect_column(c, df[c])["type"] == "numeric"]
        if len(numeric) >= 2:
            facts["layout"] = "panel"
    if x in df and y in df:
        facts["unknown_concentrations"] = bool((df[x].isna() & df[y].notna()).any())
    return facts


def open_question(key, question, allowed, proposal=None):
    """Build a public unanswered-question record.

    Args:
        key: Question id.
        question: Question definition.
        allowed: Available answers.
        proposal: Optional proposal.

    Returns:
        dict: The public record.
    """
    flag = question.get("flag")
    return {"id": key, "text": question["text"], "source": question["source"],
            "allowed": list(allowed), "proposal": proposal,
            "answer_with": f"--{flag} <answer>" if flag else f"--answer {key}=<answer>",
            "reason": "The data does not decide this answer."}


def ask_question(key, tree, answers, facts, result, node_id):
    """Resolve a question or append an open question.

    Args:
        key: Question id.
        tree: The tree.
        answers: Explicit answers.
        facts: Data facts.
        result: Mutable output.
        node_id: Current node.

    Returns:
        str or None: The answer, if resolved.
    """
    q = tree["questions"][key]
    allowed = q.get("answers", facts.get(q.get("answers_from"), []))
    explicit = answers.get(key)
    value, source, proposal = explicit, "scientist", None
    decided = None
    for rule in q.get("data_rules", []):
        if not match_rule(rule, facts):
            continue
        if "forbid" in rule and explicit in rule["forbid"]:
            raise AdviceError(rule.get("reason", f"The data rules out {key}={explicit}."),
                              "conflict", question=key, answer=explicit, fact=rule["fact"])
        if "answer" in rule:
            decided = rule["answer"]
        if "propose" in rule:
            proposal = rule["propose"]
    if key == "id" and len(facts["id.candidates"]) == 1:
        # A heuristic must not override an explicitly named valid subject column.
        proposal = facts["id.candidates"][0]["column"]
        if explicit is None:
            decided = proposal
    if explicit is not None and decided is not None and explicit != decided:
        raise AdviceError(f"The data gives {key}={decided}, not {explicit}.",
                          "conflict", question=key, answer=explicit, fact=decided)
    if explicit is None and decided is not None:
        value, source = decided, "data"
    if value is None and "default" in q:
        value, source = q["default"], "default"
    if value is None:
        result["open_questions"].append(open_question(key, q, allowed, proposal))
        return None
    if allowed and value not in allowed:
        raise AdviceError(f"Invalid {key}={value}. Choose from {', '.join(allowed)}.")
    if explicit and proposal and explicit != proposal:
        result["notes"].append(f"The data proposes {key}={proposal}; the scientist chose {explicit}.")
    if key == "design" and value == "independent" and facts["id.candidates"]:
        result["notes"].append("Subject id candidates repeat across groups: " +
                               ", ".join(c["column"] for c in facts["id.candidates"]) + ".")
    result["path"].append({"node": node_id, "question": key, "answer": value,
                           "source": source, "reason": q.get("answers", {}).get(value, q["text"]),
                           "doc_heading": tree["nodes"][node_id].get("doc_heading", result.get("doc_heading"))})
    return value


def need_flag(key, result):
    """Request a missing plot flag once.

    Args:
        key: Flag name.
        result: Mutable result.
    """
    if any(q["id"] == key for q in result["open_questions"]):
        return
    result["open_questions"].append(open_question(key, {
        "text": f"Which value should --{key} use?", "source": "scientist", "flag": key}, []))


def check_leaf_data(leaf, df, options, facts, result):
    """Check structural constraints without fitting or running a test.

    Args:
        leaf: Selected leaf.
        df: Input table.
        options: Command options.
        facts: Data facts.
        result: Mutable result.
    """
    recipe = leaf["recipe"]
    numeric_y = {"two_group_compare", "multi_group_compare", "paired_compare", "factorial_anova",
                 "feature_panel", "estimation_plot", "timecourse", "correlation", "bland_altman",
                 "dose_response", "enzyme_kinetics", "standard_curve", "distribution",
                 "forest_meta", "embedding", "manhattan", "spider_response", "waterfall"}
    if recipe in numeric_y and "y.numeric" in facts and not facts["y.numeric"]:
        raise AdviceError("The selected outcome must contain finite numeric values.", "conflict")
    if recipe in ("correlation", "bland_altman", "embedding", "manhattan", "dose_response",
                  "enzyme_kinetics", "standard_curve") and "x.numeric" in facts and not facts["x.numeric"]:
        raise AdviceError("The selected x column must contain finite numeric values.", "conflict")
    if recipe in ("survival_km", "cox_forest", "cumulative_incidence", "swimmer") and "time.numeric" in facts:
        if not facts["time.numeric"] or bool((df[options["time"]].dropna() < 0).any()):
            raise AdviceError("Time must contain finite nonnegative numbers.", "conflict")
    if recipe in ("two_group_compare", "paired_compare") and facts.get("x.n_levels", 2) != 2:
        raise AdviceError("This comparison needs exactly two groups.", "conflict")
    if recipe == "feature_panel" and options.get("test") in ("welch_t", "wilcoxon"):
        if facts.get("x.n_levels", 2) != 2:
            raise AdviceError("This feature-panel test needs exactly two groups.", "conflict")
    if recipe in ("roc_curve", "calibration") and "y.n_levels" in facts:
        if facts["y.n_levels"] != 2:
            raise AdviceError("The outcome must have exactly two levels.", "conflict")
    if recipe in ("survival_km", "cox_forest", "cumulative_incidence") and "event.codes" in facts:
        try:
            codes = [float(c) for c in facts["event.codes"]]
        except ValueError as error:
            raise AdviceError("Event codes must be numeric or logical.", "conflict") from error
        if any(not np.isfinite(c) or c < 0 or c != int(c) for c in codes):
            raise AdviceError("Event codes must be nonnegative integers.", "conflict")
        if recipe != "cumulative_incidence" and not set(codes).issubset({0, 1}):
            raise AdviceError("This recipe accepts only event codes zero and one.", "conflict")
    if recipe == "forest_meta" and not options.get("se"):
        for key in ("ci_low", "ci_high"):
            if not options.get(key):
                need_flag(key, result)
    if recipe == "feature_panel" and options.get("facet") and not options.get("y"):
        need_flag("y", result)
    if recipe == "heatmap" and options.get("test") == "welch_t" and not options.get("annotation"):
        need_flag("annotation", result)
    if options.get("annotation") and not Path(options["annotation"]).is_file():
        raise AdviceError("The annotation file does not exist.")
    if options.get("x") and options["x"] not in df and recipe not in ("roc_curve", "calibration"):
        raise AdviceError("This recipe needs a single x column.", "conflict")
    if recipe in ("roc_curve", "calibration") and options.get("x"):
        for column in split_list_option(options["x"]):
            info = column_facts(df, column, "marker")
            if not info["marker.numeric"]:
                raise AdviceError(f"Marker {column} must be finite and numeric.", "conflict")
            if recipe == "calibration" and not df[column].dropna().between(0, 1).all():
                raise AdviceError(f"Predictions in {column} must be between zero and one.", "conflict")
    if recipe in ("paired_compare", "timecourse", "sankey") or options.get("paired"):
        x, id_col = options.get("x"), options.get("id")
        if x in df and id_col in df:
            pairs = df[[x, id_col]].dropna()
            if pairs.duplicated().any():
                raise AdviceError("A subject has duplicate rows in one condition.", "conflict")
            if recipe == "timecourse":
                counts = pairs.groupby(id_col)[x].nunique()
                if bool((counts < pairs[x].nunique()).any()):
                    result["notes"].append("The repeated-measures recipe drops subjects with missing visits. "
                                           "An analysis that retains those visits needs a mixed model.")
    for flag, role in (("levels", "y" if recipe == "sankey" else "x"), ("visit_levels", "x")):
        if not options.get(flag) or f"{role}.levels" not in facts:
            continue
        levels = split_list_option(options[flag])
        if len(levels) != len(set(levels)) or set(levels) != set(facts[f"{role}.levels"]):
            raise AdviceError(f"--{flag} must list each observed level exactly once.", "conflict")
    if "x.levels" in facts:
        x = options.get("x")
        raw = list(inspect_labels(df[x]).dropna())
        if recipe in ("paired_compare", "timecourse", "sankey") or options.get("paired"):
            flag = "visit_levels" if recipe == "sankey" else "levels"
            ordered = time_level_order(raw, split_list_option(options.get(flag)))
            facts["x.plot_order"] = ordered["levels"]
        elif recipe == "proportions":
            facts["x.plot_order"] = split_list_option(options.get("levels")) or list(dict.fromkeys(raw))
        else:
            facts["x.plot_order"] = level_order(raw, split_list_option(options.get("levels")))


def finish_leaf(tree, node_id, df, options, answers, facts, result):
    """Resolve leaf requirements and produce shell-safe commands.

    Args:
        tree: The tree.
        node_id: Selected leaf node.
        df: Input table.
        options: Explicit flags.
        answers: Combined answers.
        facts: Data facts.
        result: Mutable output.
    """
    leaf = tree["nodes"][node_id]["leaf"]
    recipe = leaf["recipe"]
    result.update(recipe=recipe, doc_heading=leaf["doc_heading"])
    options = dict(options)
    for key, value in leaf.get("fixed_flags", {}).items():
        if key in options and options[key] != value:
            raise AdviceError(f"--{key} conflicts with this design.", "conflict")
        options[key] = value
    for key in leaf.get("ask_at_leaf", []):
        value = ask_question(key, tree, answers, facts, result, node_id)
        if value is not None:
            options[key] = value
    if recipe == "paired_compare":
        options.pop("paired", None)
    for key in leaf["flags"]["required"]:
        if not options.get(key):
            need_flag(key, result)
    try:
        check_recipe_options(recipe, options)
    except ValueError as error:
        raise AdviceError(str(error)) from error
    check_leaf_data(leaf, df, options, facts, result)
    test = leaf.get("test", {})
    method = options.get("test", test.get("default"))
    if method is not None and method not in test.get("allowed", []):
        raise AdviceError(f"Invalid test {method} for {recipe}. Choose from {', '.join(test.get('allowed', []))}.")
    result["test"] = {"value": method, "source": "scientist" if "test" in options else "default",
                      "checks": test.get("checks", "This recipe does not select a test from normality."),
                      "outcomes": test.get("outcomes", [])}
    for note in leaf.get("engine_notes", []):
        if note["engine"] == options.get("engine", "r"):
            result["notes"].append(note["note"])
            if method in note["forbidden"]:
                raise AdviceError(note["note"], "conflict")
    if facts.get("x.min_per_level", 3) < 3:
        result["notes"].append("A group has fewer than three rows. The recipe checks usable rows and the attainable P.")
    if result["open_questions"]:
        return
    command = [os.environ.get("PUBREADY_FIGKIT_CMD", "figkit"), "plot", "--recipe", recipe]
    for key in ("data", "engine") + tuple(sorted(set(options) - {"data", "engine"})):
        if key not in options:
            continue
        value = options[key]
        if isinstance(value, bool):
            if value:
                command.append(f"--{key}")
        else:
            # The equals form also protects values that start with a hyphen.
            if str(value).startswith("-"):
                command.append(f"--{key}={value}")
            else:
                command.extend([f"--{key}", str(value)])
    action = "plot" if "test" in options and method != "auto" else leaf.get("next_step", "plot")
    following = command.copy()
    following[1] = action
    result.update(command=command, command_text=shlex.join(command),
                  next_step={"action": action, "command": following,
                             "reason": "Show the assumption evidence to the scientist." if action == "diagnose"
                             else "The scientist can run the selected recipe."})


def advise(options, entries=(), filename=None, tree_path=None, frame=None):
    """Run the deterministic walk without any statistical computation.

    Args:
        options: Explicit plot flags.
        entries: Repeated KEY=VALUE answers.
        filename: Optional answers file.
        tree_path: Optional tree path.
        frame: Optional table for unit tests.

    Returns:
        dict: The complete public result, including its exit code.
    """
    tree, digest = load_tree(tree_path)
    result = {"advise_schema": 1, "pubready_version": os.environ.get("PUBREADY_VERSION", "3.0.0"),
              "tree_sha256": digest, "walker": "python", "status": "complete", "exit_code": 0,
              "data": None, "facts": {}, "path": [], "recipe": None, "doc_heading": None,
              "test": None, "command": None, "command_text": None, "next_step": None,
              "open_questions": [], "conflicts": [], "notes": []}
    try:
        answers = merge_answers(tree, options, entries, filename)
        options = dict(options)
        for key, value in answers.items():
            if key not in tree["questions"] or tree["questions"][key].get("flag"):
                options[key] = value
        for key, kind in PLOT_VALUE_OPTIONS:
            if key in options:
                try:
                    options[key] = kind(options[key])
                except (ValueError, TypeError) as error:
                    raise AdviceError(f"Invalid value for --{key}.") from error
        if isinstance(options.get("paired"), str):
            if options["paired"] not in ("true", "false"):
                raise AdviceError("The paired answer must be true or false.")
            options["paired"] = options["paired"] == "true"
        if options.get("engine", "r") not in ("r", "python"):
            raise AdviceError("The engine must be r or python.")
        if not options.get("data") and frame is None:
            raise AdviceError("The --data option is required.")
        df = read_tidy(options["data"], options.get("sheet")) if frame is None else frame
        result["data"] = {"file": options.get("data"), "n_rows": len(df), "n_cols": len(df.columns)}
        if df.empty:
            raise AdviceError("The data has no rows.", "conflict")
        facts = collect_facts(df, options)
        result["facts"] = facts
        # Check data-derived constraints before any missing upstream answer.
        for key in ("ref", "positive", "cause"):
            question = tree["questions"][key]
            allowed = facts.get(question["answers_from"])
            if key in answers and allowed is not None and answers[key] not in allowed:
                raise AdviceError(f"Invalid {key}={answers[key]}. Choose from {', '.join(allowed)}.")
        node_id = "start"
        while True:
            node = tree["nodes"][node_id]
            if "out_of_scope" in node:
                raise AdviceError(node["out_of_scope"], "out_of_scope")
            if "leaf" in node:
                finish_leaf(tree, node_id, df, options, answers, facts, result)
                break
            key = node["ask"]
            if key in ("outcome_type", "n_groups", "n_times", "n_factors"):
                required = "y" if key == "outcome_type" else "x"
                if not options.get(required):
                    need_flag(required, result)
                    break
                if required == "x" and facts.get("x.n_levels", 2) < 2:
                    raise AdviceError("The comparison needs at least two groups.", "conflict")
            value = ask_question(key, tree, answers, facts, result, node_id)
            if value is None:
                break
            node_id = node["branches"][value]
        if result["open_questions"]:
            result["status"] = "needs_answers"
    except AdviceError as error:
        result["status"] = error.status
        result["conflicts"].append({"reason": str(error), **error.detail})
        result["command"] = result["command_text"] = result["next_step"] = None
    result["exit_code"] = STATUS_CODES[result["status"]]
    return result


def render_result(result, output_format="json"):
    """Serialize the result with deterministic JSON or readable text.

    Args:
        result: Public result.
        output_format: json or text.

    Returns:
        str: The serialized result.
    """
    if output_format == "json":
        return json.dumps(result, ensure_ascii=False, allow_nan=False, indent=2) + "\n"
    lines = [f"figkit advise returned {result['status']}."]
    for step in result["path"]:
        lines.append(f"The answer to {step['question']} is {step['answer']} from {step['source']}.")
    if result["recipe"]:
        lines.append(f"The selected recipe is {result['recipe']}.")
    for item in result["open_questions"]:
        lines.extend([item["text"], f"Answer with {item['answer_with']}."])
        if item["allowed"]:
            lines.append("The allowed answers are " + ", ".join(item["allowed"]) + ".")
    lines.extend(item["reason"] for item in result["conflicts"])
    lines.extend(result["notes"])
    if result["command_text"]:
        lines.extend(["The plot command is:", result["command_text"],
                      result["next_step"]["reason"], shlex.join(result["next_step"]["command"])])
    return "\n".join(lines) + "\n"
