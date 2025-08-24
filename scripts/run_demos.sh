#!/usr/bin/env bash
set -euo pipefail
# 教学演示一键运行脚本
# 规则遵循：
# - 使用数组保存命令与参数
# - 双引号包裹变量
# - [[ .. ]] 条件判断

usage() {
  cat <<USAGE
用法: $(basename "$0") [--quiet|--verbose]
  --quiet    仅输出必要信息（脚本自身尽量安静）
  --verbose  额外输出脚本阶段提示
USAGE
}

main() {
  local -a cmd
  local verbose=1
  local quiet=0
  for arg in "$@"; do
    case "$arg" in
      --quiet) quiet=1; verbose=0 ;;
      --verbose) verbose=1 ;;
      -h|--help) usage; exit 0 ;;
      *) ;;
    esac
  done

  if [[ -n "${VIRTUAL_ENV:-}" ]]; then
    (( verbose )) && echo "Using venv: $VIRTUAL_ENV"
  else
    (( verbose )) && echo "Tip: 建议先创建并激活虚拟环境: python3 -m venv .venv && source .venv/bin/activate"
  fi

  cmd=(python -m examples.demo_univar)
  "${cmd[@]}"

  cmd=(python -m examples.demo_bivariate)
  "${cmd[@]}"

  cmd=(python -m examples.demo_integer_smallroot)
  "${cmd[@]}"

  cmd=(python -m examples.demo_rsa_small_e)
  "${cmd[@]}"

  cmd=(python -m examples.demo_factor_highbits)
  "${cmd[@]}"

  cmd=(python -m examples.demo_hastad_broadcast)
  "${cmd[@]}"
}

main "$@"

