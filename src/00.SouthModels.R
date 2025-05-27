modelssouth <- list(
  "Hab"=list(
    .~. + soilc,
    .~. + soilc + paspen),
  "Contrast"=list(
    .~. + road,
    .~. + road + mWell,
    .~. + road + mWell + mSoft),
  "ARU"=list(
    .~. + method),
  "Water"=list(
    .~. + pWater_KM,
    .~. + pWater_KM + pWater2_KM))
