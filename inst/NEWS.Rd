\name{msImputenews}
\title{msImpute News}
\encoding{UTF-8}

\section{Version 0.99.25}{\itemize{
    \item fix typo in \code{msImpute} man page
}}

\section{Version 0.99.24}{\itemize{
    \item Bug fix in the internal function \code{l2bary}
}}

\section{Version 0.99.23}{\itemize{
    \item \code{selectFeatures} and \code{msImpute} now use information theoretic approaches
    to find informative features for MAR/MNAR diagnosis and estimation of optimal rank, respectively.

    \item \code{lambda} in \code{msImpute} is now estimated from the data, using the bayesian interpretation of
    this shrinkage operator.

    \item \code{msImpute} can be run in three modes: "v1" is the original implementation of softImpute-als
    algorithm, "v2" is the enhanced low-rank estimation implemented in this version update, "v2-mnar"
    is adaptation of low-rank models for MNAR data. More details about methods in documentation.
    }}


\section{Version 0.99.22}{\itemize{
    \item Submitted to Bioconductor
}}


