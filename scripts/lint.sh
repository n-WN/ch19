#!/usr/bin/env bash
set -euo pipefail
# 统一 Lint/类型/测试入口（零依赖：ruff + ty + pytest）
# - 使用数组、[[ ]]、变量双引号，符合你的 shell 规范

main() {
  local -a cmd
  local do_style=0 do_types=0 do_tests=0 do_all=1

  for arg in "$@"; do
    case "$arg" in
      --style) do_style=1; do_all=0 ;;
      --types) do_types=1; do_all=0 ;;
      --tests) do_tests=1; do_all=0 ;;
      -h|--help)
        cat <<USAGE
Usage: $(basename "$0") [--style|--types|--tests]
  --style  仅运行 ruff format/check
  --types  仅运行 ty 类型检查
  --tests  仅运行 pytest
缺省运行全部步骤。
USAGE
        exit 0
        ;;
      *) ;;
    esac
  done

  if (( do_style || do_all )); then
    echo "[1/4] ruff format"
    if ! command -v ruff >/dev/null 2>&1; then
      echo "ERROR: ruff 不可用，请安装后重试" >&2; exit 1
    fi
    cmd=(ruff format .)
    "${cmd[@]}"

    echo "[2/4] ruff check (strict profile)"
    cmd=(ruff check .)
    "${cmd[@]}"
  fi

  if (( do_types || do_all )); then
    echo "[3/4] ty check (types)"
    if ! command -v uvx >/dev/null 2>&1; then
      echo "ERROR: uvx 不可用，请安装 uv/uvx 后重试" >&2; exit 1
    fi
    cmd=(uvx ty check)
    "${cmd[@]}"
  fi

  if (( do_tests || do_all )); then
    echo "[4/4] pytest"
    if ! command -v pytest >/dev/null 2>&1; then
      echo "ERROR: pytest 不可用，请安装后重试" >&2; exit 1
    fi
    cmd=(pytest)
    "${cmd[@]}"
  fi
}

main "$@"
