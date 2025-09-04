#!/usr/bin/env bash

set -euo pipefail

# ---- config ----
TEMPLATE="velmap_cum.conf"     # sample/template config containing {date}
GPS_DIR="velmap_gnss"          # folder with postseismic_gnss_YYYYMMDD.dat

# ---- checks ----
[[ -f "$TEMPLATE" ]] || { echo "ERROR: Template $TEMPLATE not found"; exit 1; }
[[ -d "$GPS_DIR"  ]] || { echo "ERROR: GPS folder $GPS_DIR not found"; exit 1; }

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

echo "Found ${#DATES[@]} date(s): ${DATES[*]}"

# ---- generate per-date configs ----
for d in "${DATES[@]}"; do
  out_cfg="velmap_cum_${d}.conf"

  # Replace {date} in the template
  sed "s/{date}/${d}/g" "$TEMPLATE" > "$out_cfg"

  echo "Wrote: $out_cfg"
done

echo "Done."

