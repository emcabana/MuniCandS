# MuniCandS

**Multivariate Tests of Uniformity, Normality and Isotropy on C and S**

## Overview

`MuniCandS` is an R package implementing Cramér-von Mises type tests for
multivariate distributions. Given an *n × p* data matrix (a sample of size *n*
in **R**^*p*), it tests whether the underlying distribution belongs to one of
the following families:

| `type` | Null hypothesis |
|--------|----------------|
| `"UC"` | Uniform on the unit hypercube [0,1]^*p* |
| `"US"` | Uniform on the hypersphere S^(*p*-1) |
| `"N"`  | Normal in **R**^*p* |
| `"I"`  | Isotropic (spherically symmetric) in **R**^*p* |
| `"E"`  | Elliptic (elliptically symmetric) in **R**^*p* |
| `"IN"` | Independent components in **R**^*p* |

The tests are based on a decomposition of a *p*-parameter Brownian sheet as
the sum of 2^*p* independent Gaussian processes, and produce two p-values
corresponding to the **m-test** and the **s-test**.

## Installation
```r
# install.packages("devtools")
devtools::install_github("emcabana/MuniCandS")
```
``` mermaid




\pagecolor{white}
\thispagestyle{empty}

\begin{center}{\Large Graphic description of the computations made by MuniCandS(Z,type)}

\begin{tikzpicture}[node distance=1.5cm]


%\node[scale=.655] (X) [ss, right of=Z,xshift=1cm] {$X$};
%\node[scale=.655] (b) [ss, right of=Z, xshift=2cm] {$\boldmath b=(b_{H^2})_{H\in\mbox{\scriptsize H\_list}}$};

\node[scale=.655] (Z1) [ss] {$Z_1$};
\node[scale=.655] (Z2) [ss, below of=Z1] {$Z_2$};
\node[scale=.655] (Zdots) [v, below of=Z2] {$\dots$};
\node[scale=.655] (Zrep) [ss, below of=Zdots] {$Z_{\mbox{\scriptsize repet}}$};

\node[scale=.655] (Z'1) [ss, below of=Zrep,yshift=-.5cm] {$Z_1^*$};
\node[scale=.655] (Z'2) [ss, below of=Z'1] {$Z_2^*$};
\node[scale=.655] (Z'dots) [v, below of=Z'2] {$\dots$};
\node[scale=.655] (Z'MC) [ss, below of=Z'dots] {$Z_{\mbox{\scriptsize MC}}^*$};

\node[scale=.655] (X1) [ss, right of=Z1,xshift=1cm] {$X_1$};
\node[scale=.655] (X2) [ss, below of=X1] {$X_2$};
\node[scale=.655] (Xdots) [v, below of=X2] {$\dots$};
\node[scale=.655] (Xrep) [ss, below of=Xdots] {$X_{\mbox{\scriptsize repet}}$};

\node[scale=.655] (b1) [ss, right of=X1,xshift=1cm] {$\boldmath b_1$};
\node[scale=.655] (b2) [ss, below of=b1] {$\boldmath b_2$};
\node[scale=.655] (bdots) [v, below of=b2] {$\dots$};
\node[scale=.655] (brep) [ss, below of=bdots] {$\boldmath b_{\mbox{\scriptsize repet}}$};

\node[scale=.655] (b'1) [ss, below of =brep,yshift=-.5cm]{$b_1^*$};
\node[scale=.655] (b'2) [ss, below of=b'1] {$b_1^*$};
\node[scale=.655] (b'dots) [v, below of=b'2] {$\dots$};
\node[scale=.655] (b'MC) [ss, below of=b'dots] {$b{\mbox{\tiny MC}}^*$};

\node[scale=.655] (X'1) [ss, below of =Xrep,yshift=-.5cm]{$X_1^*$};
\node[scale=.655] (X'2) [ss, below of=X'1] {$X_1^*$};
\node[scale=.655] (X'dots) [v, below of=X'2] {$\dots$};
\node[scale=.655] (X'MC) [ss, below of=X'dots] {$X{\mbox{\tiny MC}}^*$};


\node[scale=.655] (pvals1) [ss, right of=b'1,xshift=2.5cm]{pvals$_1$};
\node[scale=.655] (pvals2) [ss, below of=pvals1] {pvals$_2$};
\node[scale=.655] (pvalsdots) [v, below of=pvals2] {$\dots$};
\node[scale=.655] (pvalsMC) [ss, below of=pvalsdots] {pvals$_{\mbox{\tiny MC}}$};

\node[scale=.655] (pv1) [ss, right of=pvals1,xshift=1cm]{pm$_{\tiny 1}$,ps$_{\tiny 1}$};
\node[scale=.655] (pv2) [ss, below of=pv1] {pm$_{\tiny 2}$,ps$_{\tiny 2}$};
\node[scale=.655] (pvdots) [v, below of=pv2] {$\dots$};
\node[scale=.655] (pvMC) [ss, below of=pvdots] {pm$_{\mbox{\tiny MC}}$,ps$_{\mbox{\tiny MC}}$};


\node[scale=.655] (BH) [sc, right of=bdots,yshift=.7cm,xshift=-.2cm]{BH};

\node[scale=.655] (PV) [sc, right of=pvdots,yshift=.7cm,xshift=.3cm]{PV};

\node[scale=.655] (G) [sg, left of=Z'1,yshift=2.5cm,xshift=-.3cm]{G};

%\node[scale=.655] (pvals2) [ss, below of=pvals1] {$\boldmath b_2$};
%\node[scale=.655] (pvalsdots) [ss, below of=pvals2] {$\dots$};
%\node[scale=.655] (pvalsMC) [ss, below of=pvalsdots] {$\boldmath b_{\mbox{\\scriptsize MC}}$};
\node[scale=.655] (Z) [ss,below of=Z'MC,fill=blue!70!black,text=white,xshift=-.5cm,yshift=-.5cm] {arguments: $Z$, type};
\node[scale=.655] (X) [ss,below of=X'MC,yshift=-.5cm] {$X$};
\node[scale=.655] (b) [ss,below of=b'MC,yshift=-.5cm] {$b$};
\node[scale=.655] (pvals) [ss,below of=pvalsMC,yshift=-.5cm] {pvals};
\node[scale=.655] (pv) [ss,below of=pvMC,yshift=-.5cm] {pm,ps};
\node[scale=.655] (pv*) [ss,right of=pv,xshift=2cm,fill=blue!70!black,text=white] {value: pv};

\node[scale=.655] (v1) [vc, right of=b'1,xshift=.8cm]{bB2p};
\node[scale=.655] (v2) [vc, below of=v1] {bB2p};
\node[scale=.655] (vdots) [v, below of=v2] {};
\node[scale=.655] (vMC) [vc, below of=vdots] {bB2p};
\node[scale=.655] (v) [vc, below of=vMC,yshift=-.5cm] {bB2p};


\node[scale=.655] (vv) [vc, right of=pv,xshift=.3cm]{\scriptsize pP2mys};

%\draw [arrow] (Z) -- node[above] {X2b $\circ$ Z2X}(b);
%\draw[arrow] (X) --  node[above] {X2b}(b);
\draw[arrowg] (Z1) --  node[above] {\scriptsize Z2X} (X1);
\draw[arrowg] (Z2) --  node[above]  {\scriptsize  Z2X} (X2);
%\draw[arrow] (Zdots) --  node[above] {\scriptsize X2b $\circ$ Z2X} (bdots);
\draw[arrowg] (Zrep) --  node[above]  {\scriptsize  Z2X} (Xrep);
\draw[arrowg] (Z'1) --  node[above] {\scriptsize  Z2X} (X'1);
\draw[arrowg] (Z'2) --  node[above]  {\scriptsize  Z2X} (X'2);
%\draw[arrow] (Zdots) --  node[above] {\scriptsize X2b $\circ$ Z2X} (bdots);
\draw[arrowg] (Z'MC) --  node[above]  {\scriptsize  Z2X} (X'MC);
\draw[arrowg] (Z) --  node[above]  {\scriptsize  Z2X} (X);


\draw[arrow] (X1) --  node[above] {\scriptsize X2b} (b1);
\draw[arrow] (X2) --  node[above]  {\scriptsize X2b } (b2);
%\draw[arrow] (Zdots) --  node[above] {\scriptsize X2b $\circ$ Z2X} (bdots);
\draw[arrow] (Xrep) --  node[above]  {\scriptsize X2b } (brep);
\draw[arrow] (X'1) --  node[above] {\scriptsize X2b } (b'1);
\draw[arrow] (X'2) --  node[above]  {\scriptsize X2b } (b'2);
%\draw[arrow] (Zdots) --  node[above] {\scriptsize X2b $\circ$ Z2X} (bdots);
\draw[arrow] (X'MC) --  node[above]  {\scriptsize X2b } (b'MC);
\draw[arrow] (X) --  node[above]  {\scriptsize X2b } (b);


\draw[arrowg](G) -- (Z1);
\draw[arrowg](G) -- (Z2);
\draw[arrowg](G) -- (Zrep);
\draw[arrowg](G) -- (Z'1);
\draw[arrowg](G) -- (Z'2);
\draw[arrowg](G) -- (Z'MC);

\draw[arrow](b1) -- (BH);
\draw[arrow](b2) -- (BH);
%\draw[arrow](bdots) -- (BH);
\draw[arrow](brep) -- (BH);

\draw[arrow](BH) -- (v1);
\draw[arrow](BH) -- (v2);
\draw[arrow](BH) -- (vMC);
\draw[arrow](BH) -- (v);

\draw[arrow](b'1) -- (v1);
\draw[arrow](b'2) -- (v2);
\draw[arrow](b'MC) -- (vMC);
\draw[arrow](b) -- (v);

\draw[arrow](v1) -- (pvals1);
\draw[arrow](v2) -- (pvals2);
\draw[arrow](vMC) -- (pvalsMC);
\draw[arrow](v) -- (pvals);

\draw[arrow](pvals1) --  node[above]  {\scriptsize p2p}(pv1);
\draw[arrow](pvals2) -- node[above]  {\scriptsize p2p}(pv2);
\draw[arrow](pvalsMC) -- node[above]  {\scriptsize p2p}(pvMC);

\draw[arrow](pvals) -- node[above]  {\scriptsize p2p}(pv);

\draw[arrow](pv1) -- (PV);
\draw[arrow](pv2) -- (PV);
\draw[arrow](pvMC) -- (PV);
\draw[arrow](pv) -- (vv);
\draw[arrow](PV) -- (vv);
\draw[arrow](vv) -- (pv*);
\draw[arrowg](Z)-- node[left]  {\scriptsize $n,p$}(G);

\end{tikzpicture}

\end{center}
```
## Usage
```r
library(MuniCandS)

# Generate a sample from a multivariate normal distribution
set.seed(42)
X <- matrix(rnorm(200), nrow = 100, ncol = 2)

# Test normality
MuniCandS(X, type = "N")

```

## Functions

- **`MuniCandS(X, type)`** — Applies the selected test to the data matrix `X`
  and returns two p-values (m-test and s-test).

## Reference

Cabaña, A. and Cabaña, E. M. (2025). *Brownian sheet and uniformity tests on
the hypercube*. To appear in *Statistica*. arXiv:2509.06134.
<https://arxiv.org/abs/2509.06134>

## Authors

- Alejandra Cabaña — Universitat Autònoma de Barcelona, Spain
- Enrique M. Cabaña - PEDECIBA, Uruguay
