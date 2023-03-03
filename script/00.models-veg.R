mods_veg <- list(
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
        .~. + fCC2),
    "Contrast"=list(
        .~. + ROAD,
        .~. + ROAD + mWell,
        .~. + ROAD + mWell + mSoft,
        .~. + ROAD + mWell + mEnSft + mTrSft,
        .~. + ROAD + mWell + mEnSft + mTrSft + mSeism),
    "ARU"=list(
        .~. + CMETHOD),
    "Water"=list(
        .~. + pWater_KM,
        .~. + pWater_KM + pWater2_KM),
    "Space"=list(
        .~. + xPET,## climate only
        .~. + xMAT,
        .~. + xAHM,
        .~. + xFFP,
        .~. + xMAP + xFFP,
        .~. + xMAP + xFFP + xMAP:xFFP,
        .~. + xMAT + xMAP + xPET + xAHM,
        .~. + xMAT + xMAP + xPET + xAHM + xPET:xMAP + xMAT:xAHM,
        .~. + xMAT + xMAP,
        .~. + xMWMT + xMCMT,
        .~. + xAHM + xPET,
        .~.+ xY + xX + xY:xX + xPET,## linear lat-long and climate
        .~.+ xY + xX + xY:xX + xMAT,
        .~.+ xY + xX + xY:xX + xAHM,
        .~.+ xY + xX + xY:xX + xFFP,
        .~.+ xY + xX + xY:xX + xMAP + xFFP,
        .~.+ xY + xX + xY:xX + xMAP + xFFP + xMAP:xFFP,
        .~.+ xY + xX + xY:xX + xMAT + xMAP + xPET + xAHM,
        .~.+ xY + xX + xY:xX + xMAT + xMAP + xPET + xAHM + xPET:xMAP + xMAT:xAHM,
        .~.+ xY + xX + xY:xX + xMAT + xMAP,
        .~.+ xY + xX + xY:xX + xMWMT + xMCMT,
        .~.+ xY + xX + xY:xX + xAHM + xPET,
        .~.+ xY + xX + xY:xX + xY2 + xX2 + xPET,## quadratic lat-long and climate
        .~.+ xY + xX + xY:xX + xY2 + xX2 + xMAT,
        .~.+ xY + xX + xY:xX + xY2 + xX2 + xAHM,
        .~.+ xY + xX + xY:xX + xY2 + xX2 + xFFP,
        .~.+ xY + xX + xY:xX + xY2 + xX2 + xMAP + xFFP,
        .~.+ xY + xX + xY:xX + xY2 + xX2 + xMAP + xFFP + xMAP:xFFP,
        .~.+ xY + xX + xY:xX + xY2 + xX2 + xMAT + xMAP + xPET + xAHM,
        .~.+ xY + xX + xY:xX + xY2 + xX2 + xMAT + xMAP + xPET + xAHM + xPET:xMAP + xMAT:xAHM,
        .~.+ xY + xX + xY:xX + xY2 + xX2 + xMAT + xMAP,
        .~.+ xY + xX + xY:xX + xY2 + xX2 + xMWMT + xMCMT,
        .~.+ xY + xX + xY:xX + xY2 + xX2 + xAHM + xPET),
    "SSH"=list(
        .~. + SSH_KM,
        .~. + SSH05_KM),
    "HF"=list(
        .~. + THF_KM,## HF only
        .~. + Lin_KM + Nonlin_KM,
        .~. + Succ_KM + Alien_KM,
        .~. + Succ_KM + Noncult_KM + Cult_KM,
        .~. + THF_KM + THF2_KM,
        .~. + Lin_KM + Nonlin_KM + Nonlin2_KM,
        .~. + Succ_KM + Alien_KM + Succ2_KM,
        .~. + Succ_KM + Noncult_KM + Cult_KM + Succ2_KM,
        .~. + Succ_KM + Alien_KM + Alien2_KM,
        .~. + Succ_KM + Noncult_KM + Cult_KM + Noncult2_KM,
        .~. + Succ_KM + Alien_KM + Succ2_KM + Alien2_KM,
        .~. + Succ_KM + Noncult_KM + Cult_KM + Succ2_KM + Noncult2_KM),
    "Year"=list(
        .~. + YR))
