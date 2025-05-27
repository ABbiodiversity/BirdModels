modelsnorth <- list(
  "Hab"=list(
    .~. + vegc),
  "Age"=list(
    .~. + wtAge,
    .~. + wtAge + wtAge2,
    ## here dec+mix are both intercepts
    .~. + wtAge + wtAge2 + wtAge:isCon + wtAge2:isCon,
    .~. + wtAge + wtAge2 + wtAge:isUpCon + wtAge:isBogFen +
      wtAge2:isUpCon + wtAge2:isBogFen,
    .~. + wtAge + wtAge2 + wtAge:isMix + wtAge:isPine + wtAge:isWSpruce + wtAge:isBogFen +
      wtAge2:isMix + wtAge2:isPine + wtAge2:isWSpruce + wtAge2:isBogFen,
    .~. + wtAge05,
    .~. + wtAge05 + wtAge05:isCon,
    .~. + wtAge05 + wtAge05:isUpCon + wtAge05:isBogFen,
    .~. + wtAge05 + wtAge05:isMix + wtAge05:isPine + wtAge05:isWSpruce + wtAge05:isBogFen,
    .~. + wtAge05 + wtAge,
    .~. + wtAge05 + wtAge + wtAge05:isCon + wtAge:isCon,
    .~. + wtAge05 + wtAge + wtAge05:isUpCon + wtAge05:isBogFen +
      wtAge:isUpCon + wtAge:isBogFen,
    .~. + wtAge05 + wtAge + wtAge05:isMix + wtAge05:isPine + wtAge05:isWSpruce + wtAge05:isBogFen +
      wtAge:isMix + wtAge:isPine + wtAge:isWSpruce + wtAge:isBogFen),
  "CC"=list(
    .~. + fcc2),
  "Contrast"=list(
    .~. + road,
    .~. + road + mWell,
    .~. + road + mWell + mSoft,
    .~. + road + mWell + mEnSft + mTrSft,
    .~. + road + mWell + mEnSft + mTrSft + mSeism),
  "ARU"=list(
    .~. + method),
  "Water"=list(
    .~. + pWater_KM,
    .~. + pWater_KM + pWater2_KM))
