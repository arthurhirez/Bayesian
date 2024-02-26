# 1. Introduction

This work proposes a study that includes the expansion and improvement of the inverse Gompertz model - already established in statistical literature - through its uniparametric simplification. This distribution is particularly suitable for reliability and survival studies, as it models the useful life of observed phenomena. Given desirable properties, such as simplicity and efficiency, monoparametric distributions prove to be widely useful when contrasted with more complex models. The uniparametric approach reflects the flexibility of the model, promoting analyses with fewer assumptions and facilitating adaptation to a variety of data without the need to estimate multiple parameters.

Consequently, the analysis of the referred distribution is carried out through classical and Bayesian methods, with the incorporation of the paradigm of prior knowledge and uncertainties into the study, contrasting them with frequentist statistics.

# 2. The new monoparametric distribution: model construction and properties

Under the pretext of initiating the study, a brief introduction is made about the model in question: the formulation of distribution "A" is derived from the development of a cumulative distribution function, stemming from a special case of the inverse Gompertz model (a multiparametric model frequently used in demographic studies) for a single random variable.

Therefore, it is considered that a random variable X follows distribution A, characterized by the unique parameter β when its probability density function is given by:

$$ FX(x; \beta) = \exp\left(\frac{1}{\beta}\right) \left(1 - \exp\left(\frac{x}{\beta}\right)\right) $$

## 2.1 CDF Graph

![image](https://github.com/arthurhirez/Bayesiana/assets/109704516/8beb0070-78b6-48b0-93b2-eaf254fe183f)

**Figure 1:** Probability Density Function (PDF) cumulative plot created using the Seaborn library. Multiple beta values were utilized.

For \(x > 0\) and $(\beta > 0\)$. Its Probability Density Function (PDF) can consequently be found from the derivative of $F_X(x; \beta)$ concerning $x$ :

$$ f_X(x; \beta) = 1 / x^2 e^{1/\beta} (1 - e^{x/\beta}) + e^{x/\beta} \beta x $$

For $(x > 0)\$.

## 2.2 Probability Density Function (PDF) Plot
![image](https://github.com/arthurhirez/Bayesiana/assets/109704516/579a1632-96c1-45af-acb6-9888a4160238)

**Figure 2:** Plot of the probability density function.

      
