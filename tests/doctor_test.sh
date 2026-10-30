#!/usr/bin/env bash
# Tests of the Podman checks of figkit: the error for a missing Podman, for a
# Podman that cannot reach its service and for a missing image, and every result
# of `figkit doctor`. Stub `podman` and `uname` scripts on a restricted PATH
# simulate each host, so the test needs no image and changes no Podman state.
# The last case runs `figkit doctor` on the real host when its image exists.
set -uo pipefail
cd "$(dirname "$0")/.."

BASH_BIN="$(command -v bash)"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
fail=0
OUT=""
RC=0

# The commands that cli/figkit calls, other than podman and uname.
HOST_TOOLS="bash cat dirname git awk sed head tail tr grep id"

# Fill a folder with links to the host tools and a stub uname. $1 = folder.
make_bin() {
  local t
  mkdir -p "$1"
  for t in $HOST_TOOLS; do ln -s "$(command -v "$t")" "$1/$t"; done
  cat > "$1/uname" <<'EOF'
#!/usr/bin/env bash
# Stub uname: STUB_UNAME sets the kernel name (default Linux).
echo "${STUB_UNAME:-Linux}"
EOF
  chmod +x "$1/uname"
}

# Add a stub podman to a folder. $1 = folder. STUB_PODMAN sets the host state.
add_stub_podman() {
  cat > "$1/podman" <<'EOF'
#!/usr/bin/env bash
# Stub podman. STUB_PODMAN is one of: ok, rootful, no_ranges, no_service,
# no_image, mount_denied, no_machine, machine_stopped.
mode="${STUB_PODMAN:-ok}"
case "${1:-}" in --log-level=*) shift ;; esac
case "${1:-} ${2:-}" in
  "--version "*) echo "podman version 9.9.9"; exit 0 ;;
  "machine list")
    case "$mode" in
      no_machine) ;;
      machine_stopped) echo "podman-machine-default* false" ;;
      *) echo "podman-machine-default* true" ;;
    esac
    exit 0
    ;;
esac
if [ "$mode" = no_service ]; then
  # A remote client prints its own details on stdout before the error.
  [ "${1:-}" = info ] && printf 'Client:\n  OS: linux/amd64\n'
  echo "Cannot connect to Podman. Please verify your connection to the Linux system" >&2
  exit 125
fi
case "${1:-} ${2:-}" in
  "info --format")
    case "$mode" in
      rootful) echo "false 1" ;;
      no_ranges) echo "true 1" ;;
      *) echo "true 2" ;;
    esac
    ;;
  "image exists") [ "$mode" = no_image ] && exit 1; exit 0 ;;
  "image inspect") echo "sha256:stub" ;;
  "run "*)
    # STUB_LOG, when set, records the arguments of each run.
    [ -z "${STUB_LOG:-}" ] || printf '%s\n' "$*" >> "$STUB_LOG"
    if [ "$mode" = mount_denied ]; then
      echo "ls: cannot open directory '/work': Permission denied" >&2
      exit 2
    fi
    ;;
  *) echo "stub podman: unexpected call: $*" >&2; exit 99 ;;
esac
EOF
  chmod +x "$1/podman"
}

# Run cli/figkit with the given PATH; the caller sets the stub variables.
# $1 = PATH, the rest = figkit arguments. Sets OUT and RC.
run_figkit() {
  local path="$1"
  shift
  RC=0
  OUT="$(env PATH="$path" "$BASH_BIN" cli/figkit "$@" 2>&1)" || RC=$?
}

expect_rc() {  # $1 = expected exit code
  if [ "$RC" -eq "$1" ]; then
    echo "PASS: exit code $RC"
  else
    echo "FAIL: exit code $RC, expected $1"
    printf '%s\n' "$OUT" | sed 's/^/  | /'
    fail=1
  fi
}

expect_out() {  # $1 = text that the output must hold
  if printf '%s\n' "$OUT" | grep -qF -- "$1"; then
    echo "PASS: output holds \"$1\""
  else
    echo "FAIL: output lacks \"$1\""
    printf '%s\n' "$OUT" | sed 's/^/  | /'
    fail=1
  fi
}

expect_no_out() {  # $1 = text that the output must not hold
  if printf '%s\n' "$OUT" | grep -qF -- "$1"; then
    echo "FAIL: output holds \"$1\""
    printf '%s\n' "$OUT" | sed 's/^/  | /'
    fail=1
  else
    echo "PASS: output lacks \"$1\""
  fi
}

NOPODMAN="$TMP/nopodman"
STUB="$TMP/stub"
make_bin "$NOPODMAN"
make_bin "$STUB"
add_stub_podman "$STUB"
INSPECT=(inspect --data example/tumor_volume.csv)

echo "[1] no podman: a data command gives the install command, not the build hint"
run_figkit "$NOPODMAN" "${INSPECT[@]}"
expect_rc 1
expect_out "Podman is not installed"
expect_out "Install Podman with:"
expect_out "figkit doctor"
expect_no_out "Build it first"

echo "[2] no podman on macOS: the hint is brew and podman machine"
STUB_UNAME=Darwin run_figkit "$NOPODMAN" "${INSPECT[@]}"
expect_rc 1
expect_out "brew install podman && podman machine init && podman machine start"

echo "[3] no podman in Git Bash: the hint is WSL2"
STUB_UNAME=MINGW64_NT-10.0 run_figkit "$NOPODMAN" "${INSPECT[@]}"
expect_rc 1
expect_out "wsl --install -d Ubuntu"

echo "[4] no podman: figkit build stops with the install command"
run_figkit "$NOPODMAN" build
expect_rc 1
expect_out "Podman is not installed"

echo "[5] no podman: figkit version still prints the version"
run_figkit "$NOPODMAN" version
expect_rc 0
expect_out "podman: not found"

echo "[6] no podman: a command with no arguments reaches the Podman check"
run_figkit "$NOPODMAN" inspect
expect_rc 1
expect_out "Podman is not installed"
expect_no_out "unbound variable"

echo "[7] no Podman service: a data command does not ask for a build"
STUB_PODMAN=no_service run_figkit "$STUB" "${INSPECT[@]}"
expect_rc 1
expect_out "Podman cannot reach its service"
expect_out "Cannot connect to Podman"
expect_no_out "Build it first"

echo "[8] missing image: a data command asks for a build"
STUB_PODMAN=no_image run_figkit "$STUB" "${INSPECT[@]}"
expect_rc 1
expect_out "not found. Build it first: figkit build"

echo "[9] doctor, no podman"
run_figkit "$NOPODMAN" doctor
expect_rc 1
expect_out "FAIL  podman is not on the PATH"
expect_out "skip  Podman service"

echo "[10] doctor, Linux host in order"
STUB_PODMAN=ok run_figkit "$STUB" doctor
expect_rc 0
expect_out "ok    podman 9.9.9"
expect_out "ok    a test container reads the current folder"
expect_out "All checks passed."

echo "[11] doctor, rootful Podman is a warning"
STUB_PODMAN=rootful run_figkit "$STUB" doctor
expect_rc 0
expect_out "warn  Podman runs rootful"
expect_out "No check failed. Read the 1 warning(s) above."

echo "[12] doctor, no subordinate ID ranges"
STUB_PODMAN=no_ranges run_figkit "$STUB" doctor
expect_rc 1
expect_out "FAIL  no subordinate UID range"
expect_out "sudo usermod --add-subuids"

echo "[13] doctor, no Podman service"
STUB_PODMAN=no_service run_figkit "$STUB" doctor
expect_rc 1
expect_out "FAIL  Podman does not answer: Cannot connect to Podman"
expect_out "skip  rootless mode, image and test container"

echo "[14] doctor, image not built"
STUB_PODMAN=no_image run_figkit "$STUB" doctor
expect_rc 1
expect_out "is not built"
expect_out "fix: figkit build"
expect_out "skip  test container"

echo "[15] doctor, the test container cannot read the mounts"
STUB_PODMAN=mount_denied run_figkit "$STUB" doctor
expect_rc 1
expect_out "FAIL  a test container cannot read the mounted folders"
expect_out "fix: run figkit outside the sandbox"

echo "[16] doctor, macOS with no Podman machine"
STUB_UNAME=Darwin STUB_PODMAN=no_machine run_figkit "$STUB" doctor
expect_rc 1
expect_out "FAIL  no Podman machine exists"
expect_out "fix: podman machine init && podman machine start"

echo "[17] doctor, macOS with a stopped Podman machine"
STUB_UNAME=Darwin STUB_PODMAN=machine_stopped run_figkit "$STUB" doctor
expect_rc 1
expect_out "FAIL  the Podman machine is stopped"
expect_out "fix: podman machine start"

echo "[18] doctor, macOS with the folders under the home folder"
HOME="$(dirname "$(pwd -P)")" STUB_UNAME=Darwin STUB_PODMAN=ok run_figkit "$STUB" doctor
expect_rc 0
expect_out "ok    Podman machine podman-machine-default runs"
expect_out "ok    the current folder and the pubready files are in folders that the Podman machine mounts"

echo "[19] doctor, macOS with the folders outside the home folder"
HOME="$TMP" STUB_UNAME=Darwin STUB_PODMAN=ok run_figkit "$STUB" doctor
expect_rc 0
expect_out "warn  the current folder or the pubready files are outside"

echo "[20] each podman run turns off SELinux labeling"
STUB_LOG="$TMP/run.log" STUB_PODMAN=ok run_figkit "$STUB" "${INSPECT[@]}"
STUB_LOG="$TMP/run.log" STUB_PODMAN=ok run_figkit "$STUB" doctor
OUT="$(cat "$TMP/run.log" 2>/dev/null)"
expect_out "run --rm --security-opt label=disable"
if [ "$(grep -c '^run ' "$TMP/run.log")" -eq "$(grep -c 'label=disable' "$TMP/run.log")" ]; then
  echo "PASS: every recorded run has the option"
else
  echo "FAIL: a recorded run lacks --security-opt label=disable"; fail=1
fi

echo "[21] doctor on the real host"
if command -v podman >/dev/null 2>&1 && podman image exists "localhost/pubready:$(cat VERSION)" 2>/dev/null; then
  RC=0
  OUT="$(./cli/figkit doctor 2>&1)" || RC=$?
  expect_rc 0
  expect_out "ok    a test container reads the current folder"
else
  echo "SKIP: no Podman or no image on this host"
fi

echo "[22] render takes the engine from the config when no --engine is given"
printf 'engine: python\nrecipe: two_group_compare\n' > "$TMP/py_config.yaml"
printf 'engine: r\nrecipe: two_group_compare\n' > "$TMP/r_config.yaml"
STUB_LOG="$TMP/render.log" STUB_PODMAN=ok run_figkit "$STUB" render --config "$TMP/py_config.yaml"
STUB_LOG="$TMP/render.log" STUB_PODMAN=ok run_figkit "$STUB" render --config "$TMP/py_config.yaml" \
  --engine r
STUB_LOG="$TMP/render.log" STUB_PODMAN=ok run_figkit "$STUB" render --config "$TMP/r_config.yaml"
for want in "1 python3 /opt/pubready/core/python/entry/render.py" \
            "2 Rscript /opt/pubready/core/r/entry/render.R" \
            "3 Rscript /opt/pubready/core/r/entry/render.R"; do
  if sed -n "${want%% *}p" "$TMP/render.log" 2>/dev/null | grep -qF -- "${want#* }"; then
    echo "PASS: render run ${want%% *} uses ${want#* }"
  else
    echo "FAIL: render run ${want%% *} does not use ${want#* }"
    sed 's/^/  | /' "$TMP/render.log" 2>/dev/null
    fail=1
  fi
done

echo
if [ "$fail" -eq 0 ]; then
  echo "All doctor tests passed."
else
  echo "Some doctor tests failed."
fi
exit "$fail"
