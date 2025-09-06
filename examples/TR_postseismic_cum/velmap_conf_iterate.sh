#!/usr/bin/env bash
set -euo pipefail

# ---- defaults ----
TEMPLATE="velmap_cum.conf"     # template containing {date}, {insar_active}, {sboi_active}
GPS_DIR="velmap_gnss"          # has postseismic_gnss_YYYYMMDD.dat
OUT_DIR="cum_conf"             # default output dir
SBOI_MODE=0                    # 0: InSAR active, 1: SBOI active

# ---- parse flags ----
# usage: ./make_confs.sh [--sboi]
while [[ $# -gt 0 ]]; do
  case "$1" in
    --sboi) SBOI_MODE=1; shift ;;
    --template) TEMPLATE="$2"; shift 2 ;;
    --gps-dir) GPS_DIR="$2"; shift 2 ;;
    --out-dir) OUT_DIR="$2"; shift 2 ;;   # manual override if you want
    -h|--help)
      cat <<EOF
Usage: $0 [--sboi] [--template FILE] [--gps-dir DIR] [--out-dir DIR]

  --sboi          Generate SBOI configs (sets {insar_active}=0, {sboi_active}=1)
  --template      Path to template config (default: $TEMPLATE)
  --gps-dir       Directory with postseismic_gnss_*.dat (default: $GPS_DIR)
  --out-dir       Destination dir (default: $OUT_DIR; if --sboi, default becomes cum_conf_sboi)
EOF
      exit 0
      ;;
    *) echo "Unknown arg: $1" >&2; exit 2 ;;
  esac
done

# ---- adjust for SBOI mode ----
if [[ $SBOI_MODE -eq 1 && "$OUT_DIR" == "cum_conf" ]]; then
  OUT_DIR="cum_conf_sboi"
fi
INSAR_ACTIVE=$(( SBOI_MODE == 1 ? 0 : 1 ))
SBOI_ACTIVE=$(( SBOI_MODE == 1 ? 1 : 0 ))

# ---- checks ----
[[ -f "$TEMPLATE" ]] || { echo "ERROR: Template $TEMPLATE not found"; exit 1; }
[[ -d "$GPS_DIR"  ]] || { echo "ERROR: GPS folder $GPS_DIR not found"; exit 1; }
mkdir -p "$OUT_DIR"

# ---- collect dates from filenames ----
# matches: postseismic_gnss_YYYYMMDD.dat
mapfile -t DATES < <(
  find "$GPS_DIR" -maxdepth 1 -type f -name 'postseismic_gnss_*.dat' \
  | sed -E 's@.*/postseismic_gnss_([0-9]{8})\.dat@\1@g' \
  | sort -u
)

if (( ${#DATES[@]} == 0 )); then
  echo "No matching files in $GPS_DIR (postseismic_gnss_YYYYMMDD.dat). Nothing to do."
  exit 0
fi

echo "Mode: $([[ $SBOI_MODE -eq 1 ]] && echo 'SBOI' || echo 'InSAR')"
echo "INSAR_ACTIVE=$INSAR_ACTIVE  SBOI_ACTIVE=$SBOI_ACTIVE"
echo "Found ${#DATES[@]} date(s): ${DATES[*]}"

# ---- generate per-date configs ----
for d in "${DATES[@]}"; do
  out_cfg="${OUT_DIR}/velmap_cum_${d}.conf"
  # Replace placeholders: {date}, {insar_active}, {sboi_active}
  sed \
    -e "s/{date}/${d}/g" \
    -e "s/{insar_active}/${INSAR_ACTIVE}/g" \
    -e "s/{sboi_active}/${SBOI_ACTIVE}/g" \
    "$TEMPLATE" > "$out_cfg"
  echo "Wrote: $out_cfg"
done

echo "Done. Files in: $OUT_DIR/"
