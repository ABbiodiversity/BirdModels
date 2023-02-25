#' ---
#' title: "Functions"
#' author: "Peter Solymos, <solymos@ualberta.ca>"
#' date: "Nov 28, 2018"
#' output: pdf_document
#' ---
#'
#' SQL tables represent character not factor
make_char2fact <- function(x) {
    if (is.null(dim(x)))
        if (is.character(x))
            return(as.factor(x))
    for (i in seq_len(ncol(x)))
        if (is.character(x[,i]))
            x[,i] <- as.factor(x[,i])
        x
}
#' This turns (`1`, `11`) into (` 1`, `11`) [same char width]
fill_char <- function(x, fill=" ") {
    x <- as.character(x)
    n <- nchar(x)
    m <- max(n)
    sapply(seq_len(length(x)), function(i) {
        paste0(c(rep(fill, m-n[i]), x[i]), collapse="")
    })
}
#' cleans up AOU codes
normalize_sppcode <- function(x) {
    f <- function(x) {
        x <- toupper(x)
        x <- gsub(" ", "", x)
        x
    }
    if (is.factor(x)) {
        levels(x) <- f(levels(x))
        return(x)
    }
    if (is.character(x)) {
        return(f(x))
    }
    stop("x must be character or factor")
}
#' repeated calls to gsub
Gsub <- function(pattern, replacement, x, ...) {
    pattern <- as.list(pattern)
    replacement <- as.list(replacement)
    if (length(pattern) > length(replacement)) {
        j <- rep(seq_along(replacement), length(pattern))[seq_along(pattern)]
        replacement <- replacement[j]
    }
    if (length(pattern) < length(replacement)) {
        j <- rep(seq_along(pattern), length(replacement))[seq_along(replacement)]
        pattern <- pattern[j]
    }
    for (i in seq_along(pattern)) {
        x <- gsub(pattern[[i]], replacement[[i]], x, ...)
    }
    x
}
#' repeated calls to grepl, evaluating if ANY or all TRUE (NA treated as FALSE)
Grepl <- function(pattern, x, ..., ANY=TRUE) {
    if (length(pattern) == 1L)
        return(grepl(pattern, x, ...))
    z <- sapply(pattern, grepl, x=x, ...)
    z[is.na(z)] <- FALSE
    rs <- rowSums(z)
    if (ANY)
        rs > 0 else rs == ncol(z)
}
#' row standardization, avoids division by 0
row_std <- function(x) {
    rs <- rowSums(x)
    rs[rs == 0] <- 1
    #rs[is.na(rs)] <- 1
    x / rs
}
#' transforming space/climate variables
transform_clim <- function(x) {
    out <- data.frame(
        xX   = (x$X - (-113.3)) / 2.12,
        xY   = (x$Y - 54.1) / 2.43,
        xAHM = (x$AHM - 0) / 50,
        xPET = (x$PET - 0) / 800,
        xFFP = (x$FFP - 0) / 130,
        xMAP = (x$MAP - 0) / 2300,
        xMAT = (x$MAT - 0) / 6,
        xMCMT= (x$MCMT - 0) / 25,
        xMWMT= (x$MWMT - 0) / 20)
        #xASP= x$ASP
        #xSLP= log(x$SLP + 1)
        #xTRI= log(x$TRI / 5)
        #xCTI= log((x$CTI + 1) / 10)
        rownames(out) <- rownames(x)
        out$xX2 <- out$xX^2
        out$xY2 <- out$xY^2
        out
}
#' get terms from a list of formulas
get_terms <- function(mods, type=c("formula", "list"), intercept=TRUE) {
    type <- match.arg(type)
    x <- unlist(lapply(unlist(mods), function(z) as.character(z)[3]))
    #    x <- unname(substr(x, 5, nchar(x)))
    x <- gsub(". + ", "", x, fixed=TRUE)
    x <- unlist(strsplit(x, "+", fixed=TRUE))
    x <- unlist(strsplit(x, "*", fixed=TRUE))
    if (type == "list")
        x <- unlist(strsplit(x, ":", fixed=TRUE))
    x <- sub("^[[:space:]]*(.*?)[[:space:]]*$", "\\1", x, perl=TRUE)
    x <- unique(x)
    if (type == "formula") {
        x <- paste("~", paste(x, collapse=" + ", sep=""))
        if (!intercept)
            x <- paste(x, "- 1")
        x <- as.formula(x)
    }
    x
}
#' make formula term names sorted and predictable, i.e. always A:B instead of B:A
fix_names <- function(x, sep=":") {
    unlist(lapply(x, function(z) {
        paste(sort(strsplit(z, sep)[[1]]), collapse=sep)
    }))
}
#' logLik method for a failed model: make it the worst possible
'logLik.try-error' <- function (object, ...) {
    structure(-.Machine$double.xmax^(1/3), df = 1,
        nobs = 1, class = "logLik")
}
#' leave the baggage
glm_skeleton <- function(object, ..., CAICalpha=0.5) {
    ## strip parent env ref
    if (inherits(object, "try-error"))
        return(structure(as.character(object), class="try-error"))
    ## this is considered the skeleton
    out <- structure(list(
        call=object$call,
        formula=formula(object),
        coef=coef(object),
        converge=object$converge,
        logLik=as.numeric(logLik(object)),
        df=attr(logLik(object), "df"),
        nobs=nobs(object)), class="glm_skeleton")
    ## few more things that might come handy
    out$class0 <- class(object)[1L]
    out$aic <- -2*out$logLik + 2*out$df
    out$bic <- -2*out$logLik + log(out$nobs)*out$df
    out$caic <- CAICalpha * out$aic + (1-CAICalpha) * out$bic
    out
}
coef.glm_skeleton <- function(object, ...) object$coef
#' most useful stats
fstat <- function(x, level=0.95, ...) {
    .fstat <- function(x, level=0.95, ...)
        c(Mean=mean(x, ...), Median=median(x, ...),
            quantile(x, c((1-level)/2, 1 - (1-level)/2), ...))
    if (is.null(dim(x)))
        .fstat(x, level, ...) else t(apply(x, 1, .fstat, level=level, ...))
}
#' workhorse function for branching model selection
#' - `j`: bootstrap run
#' - `i`: species
#' - `mods`: model formula list
#' - `CAICalpha`: defailts to AIC
#' - `wcol`: column ID in `DAT` for weights (weights=1 if NULL)
#' - `ssh_class`: NULL or column name with land cover classes
#'   (matching column names of SSH)
#' - `ssh_fit`: model stage to be used for getting fitted values for Lc cut
#' NULL or the stage of the model to be considered for thresholding for
#'   finding 1 km$^2$ level surrounding suitable habitat.
#'   If not NULL, it expectes a stage called SSH (and predictors SSH_KM and SSH05_KM).
#' need the following objects in parent env:
#' - `DAT`: data frame with predictors
#' - `OFF`: offsets matrix
#' - `BB`: bootstrap sample ID matrix
#' - `YY`: species cross tab
#' - `SSH`: optional matrix with surrounding compositional data
.run_path1 <- function(j, i, mods, CAICalpha=1, wcol=NULL,
    ssh_class=NULL, ssh_fit=NULL)
{
    t0 <- proc.time()
    if (!is.null(ssh_class)) {
        if (!(ssh_class %in% colnames(DAT)))
            stop("ssh_class not found in DAT")
        if (is.null(ssh_fit))
            stop("specify model stage ssh_fit")
        if (!(ssh_fit %in% names(mods)))
            stop("ssh_fit not found in model stages")
        if (!("SSH" %in% names(mods)))
            stop("SSH not found in model stages")
        if (which(names(mods) == ssh_fit) >= which(names(mods) == "SSH"))
            stop("fitting stage must precede SSH stage in model list")
    }
    x <- DAT[BB[,j],]
    y <- as.numeric(YY[BB[,j], i])
    off <- if (i %in% colnames(OFF))
        OFF[BB[,j], i] else OFFmean[BB[,j]]
    w <- if (is.null(wcol))
        rep(1, nrow(x)) else x[,wcol]
    nmods <- length(mods)
    nnmods <- sapply(mods, length)
    mid <- numeric(nmods)
    bestList <- vector("list", nmods)
    caicList <- vector("list", nmods)
    names(mid) <- names(bestList) <- names(caicList) <- names(mods)
    null <- glm_skeleton(
        glm(
            formula=y ~ 1,
            data=x,
            family=poisson(),
            offset=off,
            weights=x$w,
            x=FALSE,
            y=FALSE,
            model=FALSE),
        CAICalpha=CAICalpha)
    best <- null
    res <- NULL
    for (l1 in 1:nmods) {
        if (nnmods[l1] > 0) {
            mlist <- vector("list", nnmods[l1])
            glist <- vector("list", nnmods[l1])
            for (l2 in 1:nnmods[l1]) {
                mlist[[l2]] <- glm_skeleton(try(update(object=best,
                    formula=mods[[l1]][[l2]]), silent=TRUE), CAICalpha=CAICalpha)
            }
            mcaic <- sapply(mlist, "[[", "caic")
            attr(mcaic, "StartCAIC") <- best$caic
            for (l2 in 1:length(mlist)) { # check convergence
                if (mlist[[l2]]$class != "try-error" && !mlist[[l2]]$converge)
                    mcaic[l2] <- 2*.Machine$double.xmax^(1/3)
            }
            dcaic <- mcaic - best$caic
            mmid <- which.min(dcaic)
            if (dcaic[mmid] < 0) {
                best <- mlist[[mmid]]
                mid[l1] <- mmid
                gofbest <- glist[[mmid]]
            }
            caicList[[l1]] <- mcaic
        }
        if (!is.null(ssh_class) && names(mods)[l1] == ssh_fit) {
            est <- best$coef
            Xn <- model.matrix(get_terms(mods, "formula"), x)
            colnames(Xn) <- fix_names(colnames(Xn))
            names(est) <- fix_names(names(est))
            Xn <- Xn[,names(est),drop=FALSE]
            pr <- exp(drop(Xn %*% est))

            g <- sum_by(pr, x[,ssh_class])
            l <- lorenz(g[,"x"]/g[,"by"], g[,"by"])
            s <- summary(l)
            lab <- rownames(l[l[,"x"] >= s["x[t]"],])
            x$SSH_KM <- rowSums(SSH[BB[,j],lab,drop=FALSE])
            x$SSH05_KM <- sqrt(x$SSH_KM)

            res <- list(lorenz=l, labels=lab)
        }
        bestList[[l1]] <- best
    }
    out <- list(species=i, iteration=j,
        null=null$coef,
        null_caic=null$caic,
        caic=caicList,
        coef=lapply(bestList, "[[", "coef"),
        mid=mid,
        alpha=CAICalpha,
        wcol=wcol,
        ssh=res,
        timer=proc.time()-t0)
    out
}
#' catch errors in .run_path1
run_path1 <- function(...) try(.run_path1(...))
#' bootstrap function
bfun <- function(i, SS, BLOCK=NULL) {
    set.seed(i)
    if (is.null(BLOCK))
        BLOCK <- rep(1L, length(SS))
    if (length(SS) != length(BLOCK))
        stop("lengths must be equal")
    BLOCK <- droplevels(as.factor(BLOCK))
    SS <- droplevels(as.factor(SS))
    levs <- levels(BLOCK)
    id <- seq_along(SS)
    out <- list()
    o <- sample(id)
    id1 <- id[o]
    ssyr1 <- SS[o]
    block1 <- BLOCK[o]
    for (j in levs) {
        k <- block1 == j
        ssyr0 <- ssyr1[k]
        id0 <- id1[k]
        out0 <- id0[!duplicated(ssyr0)]
        out[[j]] <- if (i == 1)
            out0 else sample(out0, length(out0), replace=TRUE)
    }
    sort(unname(unlist(out)))
}
#' Simple and fast ROC and AUC
simple_roc <- function(labels, scores){
    Labels <- labels[order(scores, decreasing=TRUE)]
    data.frame(
        TPR=cumsum(Labels)/sum(Labels),
        FPR=cumsum(!Labels)/sum(!Labels),
        Labels=Labels)
}
simple_auc <- function(ROC) {
    ROC$inv_spec <- 1-ROC$FPR
    dx <- diff(ROC$inv_spec)
    sum(dx * ROC$TPR[-1]) / sum(dx)
}
#' Making sense of model outputs
load_species <- function(path) {
    if (!file.exists(path)) {
        cat("species results file does not exist\n")
        return(NULL)
    }
    e <- new.env()
    load(path, envir=e)
    e$res
}
get_model_matrix <- function(DAT, mods) {
    X <- model.matrix(get_terms(mods, "formula"), DAT)
    colnames(X) <- fix_names(colnames(X))
    X
}
get_coef <- function(res, X, stage=NULL, na.out=TRUE) {
    OK <- !sapply(res, inherits, "try-error")
    if (any(!OK))
        warning(paste("try-error found:", sum(!OK)))
    ii <- sapply(res[OK], "[[", "iteration")
    est <- X[rep(1, length(ii)),,drop=FALSE]
    rownames(est) <- ii
    est[] <- 0
    modnams <- names(res[[ii[1]]]$coef)
    if (is.null(stage))
        stage <- length(modnams)
    if (is.character(stage)) {
        stage <- which(modnams == stage)
        if (length(stage) < 1)
            stop("stage not found")
    }
    if (stage > 0) {
        for (i in 1:length(ii)) {
            tmp <- res[[ii[i]]]$coef[[stage]]
            names(tmp) <- fix_names(names(tmp))
            sdiff <- setdiff(names(tmp), colnames(est))
            if (length(sdiff) > 0)
                stop(paste(sdiff, collapse=" "))
            est[i,match(names(tmp), colnames(est))] <- tmp
        }
    } else {
        for (i in 1:length(ii)) {
            est[i,1] <- res[[ii[i]]]$null
        }
    }
    if (any(!OK) && na.out) {
        nas <- matrix(NA, sum(!OK), ncol(est))
        rownames(nas) <- which(!OK)
        est <- rbind(est, nas)
    }
    est
}
get_caic <- function(res, stage=NULL, na.out=TRUE) {
    OK <- !sapply(res, inherits, "try-error")
    if (is.null(stage))
        stage <- length(res[[which(OK)[1]]]$coef)
    caic <- numeric(length(OK))
    caic[!OK] <- NA
    for (run in which(OK)) {
        if (stage == 0) {
            cc <- attr(res[[run]]$caic[[1]], "StartCAIC")
        } else {
            cc <- res[[run]]$caic[[stage]]
            cc <- cc[which.min(cc)]
        }
        caic[run] <- cc
    }
    if (!na.out)
        caic <- caic[OK]
    caic
}
get_summary <- function(est, show0=FALSE, ...) {
    if (!show0)
        est <- est[,colSums(abs(est), na.rm=TRUE) > 0,drop=FALSE]
    fr <- colMeans(abs(est) > 0, na.rm=TRUE)
    cf <- colMeans(est, na.rm=TRUE)
    se <- apply(est, 2, sd, na.rm=TRUE)
    z <- cf/se
    p <- 2 * pnorm(-abs(z))
    cmat <- cbind(cf, se, fr, z, p)
    colnames(cmat) <- c("Estimate", "Std. Error", "Freq.", "z value", "Pr(>|z|)")
    cmat
}
get_vcov <- function(est, show0=FALSE) {
    if (!show0)
        est <- est[,colSums(abs(est), na.rm=TRUE) > 0]
    cov(est)
}
get_confint <- function(est, level=0.95, type=c("tboot","quantile"), show0=FALSE) {
    type <- match.arg(type)
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    s <- get_summary(est, show0=show0)
    parm <- rownames(s)
    pct <- paste(format(100 * a, trim = TRUE, scientific = FALSE, digits = 3), "%", sep="")
    ci <- array(NA, dim = c(length(parm), 2), dimnames = list(parm, pct))
    if (type == "tboot") {
        fac <- qnorm(a)
        ci[] <- s[,1] + s[,2] %o% fac
    } else {
        if (!show0)
            est <- est[,colSums(abs(est)) > 0]
        cii <- t(apply(est, 2, quantile, probs=a))
        rownames(cii) <- parm
        ci[] <- cii[parm,,drop=FALSE]
    }
    return(ci)
}
#' This predicts with or without surrounding effects,
#' output is `mu` matrix (PKEY x B) that is on log scale
#' and has no offsets added to it
predict_with_SSH <- function(res, X, SSH=NULL, stage=NULL) {
    est <- suppressWarnings(get_coef(res, X, stage=stage, na.out=FALSE))
#    OK <- !sapply(res, inherits, "try-error")
#    ii <- sapply(res[OK], "[[", "iteration")
#    notNA <- which(OK)
#    est <- est[notNA,,drop=FALSE]
    c1 <- colSums(abs(est)) > 0
    if (any(c1[c("SSH_KM", "SSH05_KM")])) {
        if (is.null(SSH))
            stop("provide SSH")
        essh <- est[,c("SSH_KM", "SSH05_KM")]
        c1[c("SSH_KM", "SSH05_KM")] <- FALSE # drop SSH
        mu <- X[,c1,drop=FALSE] %*% t(est[,c1,drop=FALSE])
        mussh <- mu
        mussh[] <- 0 # put SSH effects here
        for (i in rownames(est)) {
            ssh <- res[[as.integer(i)]]$ssh
            v <- rowSums(SSH[,ssh$labels])
            mussh[,i] <- essh[1,"SSH_KM"]*v + essh[1,"SSH05_KM"]*sqrt(v)
        }
        mu <- mu + mussh # add them up
    } else {
        mu <- X[,c1,drop=FALSE] %*% t(est[,c1,drop=FALSE])
    }
    mu
}
#' Get model IDs
get_mid <- function(res) {
    OK <- !sapply(res, inherits, "try-error")
    t(sapply(res[OK], "[[", "mid"))
}
#' LinExp approximation to reduce extrapolation error in log-linear models
linexp <- function(p, beta, pmax) {
    pmax(0, 1 + p * (exp(pmax * beta) - 1) / pmax)
}
#p <- seq(0,1,0.01)
#pm <- 1/7
#b <- -1
#plot(p, exp(b*p), type="l", ylim=c(0, max(exp(b*p))))
#abline(v=pm, lty=2)
#curve(linexp(x, b, pm), add=TRUE, col=2)
#abline(h=1, lty=2)
#abline(h=linexp(1, b, pm), col=2, lty=2)
#' Trimmed mean based on quantiles
trimmed_mean <- function(x, p_range=c(0,1), ...)
    mean(x[x %[]% quantile(x, p_range[1:2], ...)], ...)
#' Uncentered correlation
#' https://stackoverflow.com/questions/23891391/uncentered-pearson-correlation
.cor_uncentered <- function(x, y) {
    sum(x*y)/(sqrt(sum(x^2)*sum(y^2)))
}
cor_uncentered <- function(x, y=NULL) {
    x <- data.matrix(x)
    y <- if (is.null(y))
        data.matrix(x) else data.matrix(y)
    if (nrow(x) != nrow(y))
        stop("rows must match")
    n <- NCOL(x)
    m <- NCOL(y)
    out <- matrix(0, n, m)
    dimnames(out) <- list(colnames(x), colnames(y))
    for (i in 1:n) {
        for (j in 1:m) {
            out[i,j] <- .cor_uncentered(x[,i], y[,j])
        }
    }
    out
}
#' Calculate Poisson CMF/PMF and deviation from 1:1
get_mass_dev <- function(x, by=0.001) {
    x <- rbind(0, x)
    Min <- pmin(x[,1], x[,2])
    Max <- pmax(x[,1], x[,2])
    Eta <- x[,"observed"]
    xout <- seq(0, 1, by=by)
    Amin <- approx(Eta, Min, xout=xout)$y
    Amax <- approx(Eta, Max, xout=xout)$y
    Amin[1] <- Amax[1] <- 0
    Amin[length(xout)] <- Amax[length(xout)] <- 1
    sum((Amax - Amin) * by)
}
get_mass <- function(y, yhat, cumulative=TRUE) {
    if (length(yhat) < length(y))
        yhat <- rep(yhat, length(y))
    yhat <- yhat[seq_len(length(y))]
    M <- max(y)
    cts <- seq(0, M)
    f <- c(sapply(cts, function(i) sum(y==i)/length(y)), 0)
    tmp <- sapply(cts, function(i) dpois(i, yhat))
    tmp <- cbind(tmp, 1-rowSums(tmp))
    fhat <- colMeans(tmp)
    if (cumulative) {
        f <- cumsum(f)
        fhat <- cumsum(fhat)
    }
    out <- cbind(observed=f, expected=fhat)
    #rownames(out) <- c(cts, paste0(">", M))
    #rownames(out) <- c(cts, paste0("(", M, ",Inf)"))
    rownames(out) <- c(cts, paste0(M+1, "+"))
    attr(out, "nobs") <- length(y)
    attr(out, "dev") <- get_mass_dev(out)
#    class(out) <- "Pmass"
    out
}
#' Monte Carlo p-value
pmc <- function(x, x0, alternative = c("two.sided", "less", "greater")) {
    alternative <- match.arg(alternative)
    pless <- sum(x >= x0, na.rm = TRUE)
    pmore <- sum(x <= x0, na.rm = TRUE)
    if (any(is.na(x0))) {
        warning("some simulated values were NA and were removed")
        nsimul <- sum(!is.na(x0))
    } else {
        nsimul <- length(x0)
    }
    p <- switch(alternative,
        two.sided = 2*pmin(pless, pmore),
        less = pless,
        greater = pmore)
    p <- pmin(1, (p + 1)/(nsimul + 1))
    p
}
#' WRSI function
wrsi <-
function(y, x)
{
    Y <- as.integer(ifelse(y > 0, 1L, 0L))
    X <- data.matrix(x)
    n <- length(Y)
    if (nrow(X) != n)
        stop("dimension mismatch: X and Y")
    K <- ncol(X)
    if (is.null(colnames(X)))
        colnames(X) <- paste0("V", seq_len(K))
    ## relative suitability index
    ## # of available units of type k / total # of available units (any type)
    Pavail <- colSums(X) / sum(X)
    ## # of used units of type k / total # of used units (any type)
    Xu <- X * Y
    ## sum(Xu) = sum(Y) except when rowsum != 1
    Pused <- colSums(Xu) / sum(Xu)
    ## crude weighted p-occ
    Pw <- colSums(Xu) / colSums(X)
    ## Weighted Relative Suitability Index
    WRSI <- Pused / Pavail
    #Var <- (1/colSums(Xu)) - (1/sum(Xu)) + (1/colSums(X)) - (1/sum(X))
    res <- data.frame(
        WRSI=WRSI,
        zWRSI=log(WRSI),
        rWRSI=(exp(2 * log(WRSI)) - 1)/(1 + exp(2 * log(WRSI))),
        Pused=Pused,
        Pavail=Pavail,
        Pw=Pw,
        u=colSums(Xu),
        a=colSums(X))
    rownames(res) <- colnames(X)
    class(res) <- c("wrsi", "data.frame")
    res
}

ssh_freq <- function(res) {
    LAB <- sort(rownames(res[[1]]$ssh$lorenz)[-1L])
    lab <- list()
    for (i in seq_len(length(res))) {
        if (!inherits(res[[i]], "try-error"))
            lab[[length(lab)+1L]] <- res[[i]]$ssh$labels
    }
    labv <- factor(unlist(lab), LAB)
    o <- data.frame(table(Label=labv))
    o$Prop <- o$Freq / length(lab)
    o[order(o$Prop, decreasing = TRUE),]
}


get_coef_north <- function(est, subset=1, ...) {
    estn <- est[subset,]

    mu <- drop(Xage %*% estn[colnames(Xage)])
    lam1 <- exp(mu)
    lam1 <- lam1[!grepl("9", names(lam1))]
    lamCC <- lam1[grepl("CC", names(lam1))]

    #    MOD <- c("ROAD", "mWell", "mSoft", "mEnSft", "mTrSft", "mSeism")
    MOD <- c("mWell", "mSoft", "mEnSft", "mTrSft", "mSeism")
    Z <- exp(estn[MOD])
    isSoft <- estn["mSoft"] != 0 & estn["mEnSft"] == 0
    #isSoft2 <- get_mid(resn)[,"Contrast"] == 3
    if (isSoft) {
        estn["mEnSft"] <- estn["mSoft"]
        estn["mTrSft"] <- estn["mSoft"]
        estn["mSeism"] <- estn["mSoft"]
    }
    # pm <- c("ROAD"=1, "mWell"=0.2, "mSoft"=0.2,
    #         "mEnSft"=0.2, "mTrSft"=0.2, "mSeism"=0.05)
    # pm <- c("ROAD"=1, "mWell"=0.1, "mSoft"=0.1,
    #         "mEnSft"=0.1, "mTrSft"=0.1, "mSeism"=0.025)
    # for (i in MOD)
    #     Z[i] <- linexp(1, estn[i], pm[i])

    HFc <- c("Crop", "Industrial", "Mine", "RoughP", "Rural", "TameP", "Urban")

    FUN <- mean

    # not in HFc and not forestry!
    Xn2 <- Xn[en$DAT$mWell > 0 & en$DAT$vegc %ni% HFc & en$DAT$fCC2 == 0, colnames(Xage)]
    nWell <- nrow(Xn2)
    lamWell <- exp(Xn2 %*% estn[colnames(Xn2)])
    estWell <- FUN(lamWell * Z["mWell"])

    Xn2 <- Xn[en$DAT$mEnSft > 0 & en$DAT$vegc %ni% HFc & en$DAT$fCC2 == 0, colnames(Xage)]
    nEnSoft <- nrow(Xn2)
    lamEnSft <- exp(Xn2 %*% estn[colnames(Xn2)])
    estEnSft <- FUN(lamEnSft * Z["mEnSft"])

    ## TrSft incorporates ROAD effect as well?
    Xn2 <- Xn[en$DAT$mTrSft > 0 & en$DAT$vegc %ni% HFc & en$DAT$fCC2 == 0, colnames(Xage)]
    nTrSoft <- nrow(Xn2)
    lamTrSft <- exp(Xn2 %*% estn[colnames(Xn2)])
    #estTrSft <- quantile(lamTrSft * Z[,"mTrSft"] * Z[,"ROAD"], c(0.5, 0.05, 0.95))
    estTrSft <- FUN(lamTrSft * Z["mTrSft"])

    Xn2 <- Xn[en$DAT$mSeism > 0 & en$DAT$vegc %ni% HFc & en$DAT$fCC2 == 0, colnames(Xage)]
    nSeism <- nrow(Xn2)
    lamSeism <- exp(Xn2 %*% estn[colnames(Xn2)])
    estSeism <- FUN(lamSeism * Z["mSeism"])


    ## capping values
    ## early seral: young forest (R) and open/shrub
    EARLY <- c(names(lam1)[endsWith(names(lam1), "R")],
               "GrassHerb", "Shrub", "GraminoidFen", "Marsh")
    ESMAX <- max(lam1[EARLY])
    ## non vegetated HF: this should apply to wellsites besides early seral
    NONVEG <-  c("Industrial", "Rural", "Urban")
    NVMAX <- max(lam1[NONVEG], ESMAX)

    ## UrbInd
    #Xn2 <- Xn[en$DAT$vegc %ni% c("Industrial", "Mine", "Rural", "Urban"), colnames(Xage)]
    #nUI <- nrow(Xn2)
    #lamUI <- apply(exp(Xn2 %*% t(estn[,colnames(Xn2),drop=FALSE])), 2, median)
    #estUI <- quantile((nUI * lamUI + nWell * lamWell) / (nUI + nWell),
    #    c(0.5, 0.05, 0.95))

    lam1 <- c(lam1,
              Wellsites=min(NVMAX, estWell), # nonveg HF + early seral
              EnSeismic=min(ESMAX, estSeism),# early seral
              EnSoftLin=min(ESMAX, estEnSft),# early seral
              TrSoftLin=min(ESMAX, estTrSft),# early seral
              HardLin=0,
              Water=0,
              Bare=0,
              SnowIce=0,
              MineV=unname(lam1["Mine"]))
    lam1["Mine"] <- 0
    names(lam1) <- gsub("Spruce", "WhiteSpruce", names(lam1))
    names(lam1) <- gsub("Decid", "Deciduous", names(lam1))
    lam1
}

get_coef_south <- function(est, subset=1, ...) {
    ests <- est[subset,]
    LCC <- c("soilcBlowout", "soilcClaySub", "soilcCrop",
             "soilcIndustrial", "soilcMine", "soilcOther", "soilcRapidDrain",
             "soilcRoughP", "soilcRural", "soilcSandyLoam", "soilcTameP",
             "soilcThinBreak", "soilcUrban")
    LCC2 <- c("soilc2Crop", "soilc2OtherBlowThinRapid",
              "soilc2RoughP", "soilc2Rural", "soilc2TameP", "soilc2UrbInd")
    LCC1 <- c("soilc1Crop", "soilc1RoughP", "soilc1TameP", "soilc1UrbIndRur")
    if (all(ests[LCC] == 0)) {
        ## soilc2 in use
        if (all(ests[LCC1] == 0)) {
            #ests[c("Loamy", "ClaySub", "SandyLoam")] <- ests["ClaySubLoamSand"] # intercept
            ests[c("soilcBlowout", "soilcOther", "soilcRapidDrain",
                   "soilcThinBreak")] <- ests["soilc2OtherBlowThinRapid"]
            ests[c("soilcIndustrial", "soilcMine", "soilcUrban")] <- ests["soilc2UrbInd"]
            ests[c("soilcCrop", "soilcRoughP", "soilcRural", "soilcTameP")] <-
                ests[c("soilc2Crop", "soilc2RoughP", "soilc2Rural", "soilc2TameP")]
        } else {
            ## soilc1 in use
            #            ests[c("soilcBlowout", "soilcOther", "soilcRapidDrain", "soilcThinBreak",
            #                "soilcLoamy", "soilcClaySub", "soilcSandyLoam")] <- ests["soilc1SoilNative"]
            ests[c("soilcIndustrial", "soilcMine", "soilcUrban",
                   "soilcRural")] <- ests["soilc1UrbIndRur"]
            ests[c("soilcCrop", "soilcRoughP", "soulcTameP")] <-
                ests[c("soilc1Crop", "soilc1RoughP", "soulc1TameP")]
        }
    }
    muLCC <- c(ests[1], ests[1]+ests[LCC])
    names(muLCC) <- levels(es$DAT$soilc)
    lamLCC0 <- exp(muLCC) # nontreed
    o <- c("Loamy", "Blowout", "ClaySub", "RapidDrain", "SandyLoam", "ThinBreak", "Other",
           "RoughP", "TameP","Crop","Urban", "Rural", "Industrial", "Mine")
    lamLCC0 <- lamLCC0[o]

    #    MOD <- c("ROAD", "mWell", "mSoft")
    MOD <- c("mWell", "mSoft")
    Z <- exp(ests[MOD])
    # pm <- c("ROAD"=1, "mWell"=0.2, "mSoft"=0.2)
    # for (i in MOD)
    #     Z[i] <- linexp(1, ests[i], pm[i])
    lamMOD <- Z
    #    names(lamMOD) <- c("Road", "Well", "Soft")
    names(lamMOD) <- c("Well", "Soft")

    HFc <- c("RoughP", "TameP","Crop","Urban", "Rural", "Industrial", "Mine")

    FUN <- mean

    Xs2 <- Xs[es$DAT$mSoft > 0 & es$DAT$soilc %ni% HFc, colnames(Xs)]
    nSoft <- nrow(Xs2)
    lamSoft <- exp(Xs2 %*% ests[colnames(Xs2)])
    estSoft <- FUN(lamSoft * Z["mSoft"])

    Xs2 <- Xs[es$DAT$mWell > 0 & es$DAT$soilc %ni% HFc, colnames(Xs)]
    nWell <- nrow(Xs2)
    lamWell <- exp(Xs2 %*% ests[colnames(Xs2)])
    estWell <- FUN(lamWell * Z["mWell"])

    #    Xs2 <- Xs[es$DAT$ROAD > 0 & es$DAT$soilc %ni% HFc, colnames(Xs)]
    #    nROAD <- nrow(Xs2)
    #    lamROAD <- exp(Xs2 %*% ests[colnames(Xs2)])
    #    estROAD <- FUN(lamROAD * Z["ROAD"])

    ## capping values
    ## early seral: young forest (R) and open/shrub
    EARLY <- c("Loamy", "Blowout", "ClaySub", "RapidDrain", "SandyLoam", "ThinBreak", "Other")
    ESMAX <- max(lamLCC0[EARLY])
    ## non vegetated HF: this should apply to wellsites besides early seral
    NONVEG <-  c("Industrial", "Rural", "Urban")
    NVMAX <- max(lamLCC0[NONVEG], ESMAX)

    lam0 <- c(lamLCC0,
              Wellsites=min(NVMAX, estWell),
              EnSeismic=min(ESMAX, estSoft),
              EnSoftLin=min(ESMAX, estSoft),
              TrSoftLin=min(ESMAX, estSoft),
              HardLin=0, # survace is unsuitable as hell
              #HFor=0,
              MineV=unname(lamLCC0["Mine"]),
              Water=0)
    lam0["Mine"] <- 0
    attr(lam0, "pAspen") <- ests["pAspen"]
    lam0
}

ALLBIRDSPP <- list(
    north = c(AlderFlycatcher = "ALFL", AmericanCrow = "AMCR", AmericanGoldfinch = "AMGO",
              AmericanRedstart = "AMRE", AmericanRobin = "AMRO", AmericanThreetoedWoodpecker = "ATTW",
              BaltimoreOriole = "BAOR", BarnSwallow = "BARS", BarredOwl = "BADO", BlackAndWhiteWarbler = "BAWW",
              BlackbilledMagpie = "BBMA", BaybreastedWarbler = "BBWA", BlackbackedWoodpecker = "BBWO",
              BlackcappedChickadee = "BCCH", BrownheadedCowbird = "BHCO", BlueheadedVireo = "BHVI",
              BlackburnianWarbler = "BLBW", BlueJay = "BLJA", BlackpollWarbler = "BLPW",
              BorealChickadee = "BOCH", BrewersBlackbird = "BRBL", BrownCreeper = "BRCR",
              BlackthroatedGreenWarbler = "BTNW", CanadaWarbler = "CAWA", ClarksNutcracker = "CLNU", ClaycoloredSparrow = "CCSP",
              CedarWaxwing = "CEDW", ChippingSparrow = "CHSP", CapeMayWarbler = "CMWA", CommonLoon="COLO",
              CommonNighthawk = "CONI", ConnecticutWarbler = "CONW", CommonRaven = "CORA",
              CommonYellowthroat = "COYE", ChestnutsidedWarbler = "CSWA", DarkeyedJunco = "DEJU",
              DownyWoodpecker = "DOWO", EasternKingbird = "EAKI", EasternPhoebe = "EAPH",
              EuropeanStarling = "EUST", EveningGrosbeak = "EVGR", FoxSparrow = "FOSP",
              GoldencrownedKinglet = "GCKI", GrayJay = "GRAJ", GrayCatbird = "GRCA",
              GreaterYellowlegs = "GRYE", HairyWoodpecker = "HAWO", HermitThrush = "HETH", HornedGrebe = "HOGR",
              HornedLark = "HOLA", HouseSparrow = "HOSP", HouseWren = "HOWR",
              Killdeer = "KILL", LeContesSparrow = "LCSP", LeastFlycatcher = "LEFL",
              LesserYellowlegs = "LEYE", LincolnsSparrow = "LISP", MagnoliaWarbler = "MAWA",
              MarshWren = "MAWR", MountainChickadee = "MOCH", MourningDove = "MODO",
              MourningWarbler = "MOWA", NelsonsSparrow = "NESP", NorthernFlicker = "NOFL", NorthernPintail = "NOPI",
              NorthernSawWhetOwl = "NSWO", NorthernWaterthrush = "NOWA", OrangecrownedWarbler = "OCWA",
              OlivesidedFlycatcher = "OSFL", Ovenbird = "OVEN", PalmWarbler = "PAWA",
              PhiladelphiaVireo = "PHVI", PineSiskin = "PISI", PileatedWoodpecker = "PIWO",
              PurpleFinch = "PUFI", RedTailedHawk = "RTHA", RosebreastedGrosbeak = "RBGR", RedbreastedNuthatch = "RBNU",
              RubycrownedKinglet = "RCKI", RedeyedVireo = "REVI", RockPigeon = "ROPI",
              RustyBlackbird = "RUBL", RuffedGrouse = "RUGR", RedwingedBlackbird = "RWBL",
              SavannahSparrow = "SAVS", SedgeWren = "SEWR", Sora = "SORA",
              SolitarySandpiper = "SOSA", SongSparrow = "SOSP", SpottedSandpiper = "SPSA",
              SwampSparrow = "SWSP", SwainsonsThrush = "SWTH", TennesseeWarbler = "TEWA",
              TreeSwallow = "TRES", VariedThrush = "VATH", Veery = "VEER",
              VesperSparrow = "VESP", WarblingVireo = "WAVI", WhitebreastedNuthatch = "WBNU",
              WhitecrownedSparrow = "WCSP", WesternMeadowlark = "WEME", WesternTanager = "WETA",
              WesternWoodPewee = "WEWP", WilsonsSnipe = "WISN", WilsonsWarbler = "WIWA",
              WinterWren = "WIWR", WhitethroatedSparrow = "WTSP", WhitewingedCrossbill = "WWCR",
              YellowbelliedFlycatcher = "YBFL", YellowbelliedSapsucker = "YBSA", YellowRail = "YERA",
              YellowWarbler = "YEWA", YellowheadedBlackbird = "YHBL", YellowrumpedWarbler = "YRWA"),
    south = c(AlderFlycatcher = "ALFL", AmericanCrow = "AMCR", AmericanGoldfinch = "AMGO",
              AmericanRobin = "AMRO", BairdsSparrow = "BAIS", BankSwallow = "BANS",
              BaltimoreOriole = "BAOR", BarnSwallow = "BARS", BlackbilledMagpie = "BBMA",
              BrownheadedCowbird = "BHCO", Bobolink = "BOBO", BrewersBlackbird = "BRBL",
              BrewersSparrow = "BRSP", BrownThrasher = "BRTH", ChestnutcollaredLongspur = "CCLO",
              ClaycoloredSparrow = "CCSP", CedarWaxwing = "CEDW", ChippingSparrow = "CHSP",
              CommonGrackle = "COGR", CommonNighthawk = "CONI", CommonRaven = "CORA", CommonYellowthroat = "COYE",
              DownyWoodpecker = "DOWO", EasternKingbird = "EAKI", EasternPhoebe = "EAPH",
              EuropeanStarling = "EUST", GrayPartridge = "GRAP", GrasshopperSparrow = "GRSP",
              HornedLark = "HOLA", HouseSparrow = "HOSP", HouseWren = "HOWR",
              Killdeer = "KILL", LarkBunting = "LARB", LarkSparrow = "LASP",
              LongbilledCurlew = "LBCU", LeContesSparrow = "LCSP", LeastFlycatcher = "LEFL",
              LincolnsSparrow = "LISP", MarbledGodwit = "MAGO", McCownsLongspur = "MCLO",
              MourningDove = "MODO", NelsonsSparrow = "NESP", NorthernFlicker = "NOFL", NorthernPintail = "NOPI",
              RingneckedPheasant = "RNEP", RockPigeon = "ROPI", RedTailedHawk = "RTHA", RedwingedBlackbird = "RWBL",
              SavannahSparrow = "SAVS", Sora = "SORA", SongSparrow = "SOSP",
              SpraguesPipit = "SPPI", SpottedSandpiper = "SPSA", SharptailedGrouse = "STGR",
              TreeSwallow = "TRES", UplandSandpiper = "UPSA", VesperSparrow = "VESP",
              WarblingVireo = "WAVI", WesternKingbird = "WEKI", WesternMeadowlark = "WEME",
              WesternWoodPewee = "WEWP", Willet = "WILL", WilsonsPhalarope = "WIPH",
              WilsonsSnipe = "WISN", WhitethroatedSparrow = "WTSP", YellowWarbler = "YEWA",
              YellowheadedBlackbird = "YHBL"))
