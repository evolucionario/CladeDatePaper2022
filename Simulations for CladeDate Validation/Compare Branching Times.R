###############################
### Compare Branching Times ###
###############################

# Function modified from function comparePhylo() in the 'ape' package (Paradis & Schielp 2021)
# A simplified version comparePhylo() focused on comparing branching times.

# Paradis, E., & Schliep, K. (2019). ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. Bioinformatics, 35(3), 526-528.


compareBTs <- function (x, y, force.rooted = FALSE) 
{
    tree1 <- deparse(substitute(x))
    tree2 <- deparse(substitute(y))
    res <- list()
    msg <- paste("=> Comparing", tree1, "with", tree2)
    res$messages <- msg
    n1 <- Ntip(x)
    n2 <- Ntip(y)
    tmp <- if (n1 == n2) 
        paste("Both trees have the same number of tips:", n1)
    else paste("Trees have different numbers of tips:", n1, "and", 
        n2)
    msg <- c(msg, tmp)
    tips1 <- x$tip.label
    tips2 <- y$tip.label
    tips12 <- match(tips1, tips2)
    tips21 <- match(tips2, tips1)
    tmp <- is.na(tips12)
    if (any(tmp)) 
        msg <- c(msg, paste("Tips in", tree1, "not in", tree2, 
            ":", paste(tips1[tmp], collapse = ", ")))
    tmp2 <- is.na(tips21)
    if (any(tmp2)) 
        msg <- c(msg, paste("Tips in", tree2, "not in", tree1, 
            ":", paste(tips2[tmp2], collapse = ", ")))
    sameTips <- FALSE
    if (!sum(tmp, tmp2)) {
        msg <- c(msg, "Both trees have the same tip labels")
        sameTips <- TRUE
    }
    m1 <- Nnode(x)
    m2 <- Nnode(y)
    tmp <- if (m1 == m2) 
        paste("Both trees have the same number of nodes:", m1)
    else paste("Trees have different numbers of nodes:", m1, 
        "and", m2)
    msg <- c(msg, tmp)
    rooted1 <- is.rooted(x)
    rooted2 <- is.rooted(y)
    tmp <- if (rooted1) {
        if (rooted2) 
            "Both trees are rooted"
        else paste(tree1, "is rooted,", tree2, "is unrooted")
    }
    else {
        if (rooted2) 
            paste(tree1, "is unrooted,", tree2, "is rooted")
        else "Both trees are unrooted"
    }
    msg <- c(msg, tmp)
    ultra1 <- ultra2 <- FALSE
    if (!is.null(x$edge.length)) 
        ultra1 <- is.ultrametric(x)
    if (!is.null(y$edge.length)) 
        ultra2 <- is.ultrametric(y)
    tmp <- if (ultra1) {
        if (ultra2) 
            "Both trees are ultrametric"
        else paste(tree1, "is ultrametric,", tree2, "is not")
    }
    else {
        if (ultra2) 
            paste(tree1, "is not ultrametric,", tree2, "is ultrametric")
        else "Both trees are not ultrametric"
    }
    msg <- c(msg, tmp)
    if (rooted1 && rooted2) {
        key1 <- makeNodeLabel(x, "md5sum")$node.label
        key2 <- makeNodeLabel(y, "md5sum")$node.label
        mk12 <- match(key1, key2)
        mk21 <- match(key2, key1)
        if (any(tmp <- is.na(mk12))) {
            nk <- sum(tmp)
            msg <- c(msg, paste(nk, if (nk == 1) "clade" else "clades", 
                "in", tree1, "not in", tree2))
        }
        if (any(tmp <- is.na(mk21))) {
            nk <- sum(tmp)
            msg <- c(msg, paste(nk, if (nk == 1) "clade" else "clades", 
                "in", tree2, "not in", tree1))
        }
        nodes1 <- which(!is.na(mk12))
        nodes2 <- mk12[!is.na(mk12)]
        if (ultra1 && ultra2) {
            bt1 <- branching.times(x)
            bt2 <- branching.times(y)
            BT <- data.frame(bt1[nodes1], bt2[nodes2])
            names(BT) <- c(tree1, tree2)
            res$BT <- BT
            msg <- c(msg, "Branching times of clades in common between both trees: see ..$BT\n(node number in parentheses)")
        }
        if (!is.null(nl1 <- x$node.label) && !is.null(nl2 <- y$node.label)) {
            NODES <- data.frame(paste0(nl1[nodes1], " (", nodes1 + 
                n1, ")"), paste0(nl2[nodes2], " (", nodes2 + 
                n2, ")"))
            names(NODES) <- c(tree1, tree2)
            res$NODES <- NODES
            msg <- c(msg, "Node labels of clades in common between both trees: see ..$NODES\n(node number in parentheses)")
        }
    }
    if (!force.rooted && !rooted1 && !rooted2 && sameTips && 
        m1 == m2) {
        TR <- .compressTipLabel(c(x, y))
        bs <- bitsplits(TR)
        common.splits <- which(bs$freq == 2L)
        ncs <- length(common.splits)
        tmp <- if (ncs) 
            paste(ncs, if (ncs == 1) 
                "split"
            else "splits", "in common")
        else "No split in common"
        msg <- c(msg, tmp)
    }
    res$messages <- paste0(msg, ".")
    class(res) <- "comparePhylo"
    res
}
