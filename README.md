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

# 2.3 Properties
## 2.3.1 Unimodality
From the graphical analysis of the probability density functions (see Figure 1), the unimodal behavior of the PDF is evident, which is explored in this subsection. It is observed that the data is concentrated in a single peak, facilitating the perception of its central distribution. The mode is calculated through the derivative of the PDF using the following equation:

$$ \frac{d}{dx} f_X(x; \beta) = 0 \quad \text{exp} \left( \frac{x}{\beta} \right) - 2 \frac{x}{\beta} - 1 = 0 $$

## 2.3.2 Reliability
Reliability is defined as the probability of a system continuing to operate without failure until a certain point in time. The reliability function can be obtained from the complement of the cumulative distribution function (CDF). Thus, the reliability function \( R_X(x; \beta) \) is expressed as:

$$ R_X(x; \beta) = 1 - \text{exp} \left( \frac{1}{\beta} \left( 1 - \text{exp} \left( \frac{x}{\beta} \right) \right) \right) $$

![image](https://github.com/arthurhirez/Bayesiana/assets/109704516/cc0148f4-6544-4cb2-92b7-0a8db471685d)

**Figure 3:** Plot of the reliability function

## 2.3.3 Hazard Function
The hazard rate of a given system is defined as the probability of "instantaneous failure" given that it has "survived" until a given moment. The hazard function can be obtained by dividing the PDF by the complement of the CDF:

$$ h(t) = \frac{f(t)}{1 - F(t)} $$

$$ h_X(x; \beta) = \frac{\exp \left( \frac{x}{\beta} \right) x^2}{\exp \left( \frac{1}{\beta} \right) \exp \left( \frac{x}{\beta} \right) - 1} $$

![image](https://github.com/arthurhirez/Bayesiana/assets/109704516/948b0a03-a885-433f-b355-b3f9f2a8d318)

**Figure 4:** Plot of the hazard function

It is noteworthy for its unconventional shape, with an initially high failure rate followed by a decrease that stabilizes. It suggests that for a given system, after overcoming the critical period, it becomes progressively more reliable.


# 2.3.4 Quantile Function
The quantile function of the distribution, used for generating samples in simulation studies, can essentially be obtained through the inverse of the CDF. In this study, the quantile function is derived using numerical methods. Its analytical expression is given by:

$$ x_q = \frac{\ln(1 - \beta \ln(q))}{\beta} $$

## 3 Simulation Study using the Classical Method
To evaluate the effectiveness of the Maximum Likelihood Estimation method for the parameter \( \beta \), a simulation study was conducted on samples generated in multiple sizes. Each step of the experiment is described in this section.

### 3.1 Data Generation
To standardize the experiments for more reliable results, a file containing a thousand replicas of uniform data was generated for each size explored in the simulation. Using the inverse transformation method, i.e., through the quantile function of the distribution, random samples from the distribution under analysis can be obtained by calculating \( X_i = F^{-1}(U_i) \) for the previously generated data.

### 3.2 Maximum Likelihood Estimator
Given the previously obtained samples, the likelihood function is the product of the density functions for each data point, depending on the distribution parameter. Therefore:

$$ L(\beta; x_1, x_2, ..., x_n) = \prod_{i=1}^{n} f_X(x_i; \beta) $$

Subsequently, the logarithmic function is applied to the obtained expression to simplify the calculations by transforming products into sums without altering the behavior of the function, i.e., the global maximizer of the likelihood remains unchanged.

The goal is to find the value of \( \beta \) that maximizes the log-likelihood, obtained from the first derivative of the equation set to zero:

$$ LL = \frac{1}{\beta} - \sum_{i=1}^{n} \frac{x_i}{\beta \exp(\beta x_i)} - \sum_{i=1}^{n} (1 - 2 \ln(x_i)) $$

The solution to the equation is approximated using numerical methods due to the impossibility of an analytical solution. In this case, the "fsolve" function from the scipy library leverages nonlinear optimization concepts to approximate the solution of the equation without its explicit derivative.

The objective is to estimate the bias and variance metrics of the estimator obtained through the maximum likelihood method. Specifically, the RAB, mean, MSE, RMSE, and coverage probability are calculated for each replication.


**Table 1:** Results of Simulation for β = 0.125
| Size | Mean     | Var      | Bias     | MSE      | RMSE     | PC    |
|------|----------|----------|----------|----------|----------|-------|
| 25   | 0.214629 | 0.044276 | 0.089629 | 0.052309 | 0.228712 | 0.933 |
| 50   | 0.173871 | 0.019995 | 0.048871 | 0.022384 | 0.149612 | 0.935 |
| 100  | 0.145898 | 0.007447 | 0.020898 | 0.007884 | 0.088791 | 0.945 |
| 200  | 0.136500 | 0.003540 | 0.011500 | 0.003672 | 0.060599 | 0.947 |
| 400  | 0.131363 | 0.001713 | 0.006363 | 0.001754 | 0.041880 | 0.948 |


**Table 2:** Results of Simulation for β = 0.6
| Size | Mean     | Var      | Bias     | MSE      | RMSE     | PC    |
|------|----------|----------|----------|----------|----------|-------|
| 25   | 0.698618 | 0.074206 | 0.098618 | 0.083932 | 0.289710 | 0.932 |
| 50   | 0.653746 | 0.035250 | 0.053746 | 0.038139 | 0.195292 | 0.940 |
| 100  | 0.622565 | 0.013778 | 0.022565 | 0.014287 | 0.119529 | 0.951 |
| 200  | 0.612607 | 0.006740 | 0.012607 | 0.006899 | 0.083057 | 0.949 |
| 400  | 0.607330 | 0.003317 | 0.007330 | 0.003370 | 0.058056 | 0.946 |


**Table 3:** Results of Simulation for β = 1
| Size | Mean     | Var      | Bias     | MSE      | RMSE     | PC    |
|------|----------|----------|----------|----------|----------|-------|
| 25   | 1.106558 | 0.100101 | 0.106558 | 0.111455 | 0.333849 | 0.932 |
| 50   | 1.058147 | 0.048430 | 0.058147 | 0.051811 | 0.227621 | 0.941 |
| 100  | 1.024350 | 0.019263 | 0.024350 | 0.019856 | 0.140910 | 0.954 |
| 200  | 1.013706 | 0.009520 | 0.013706 | 0.009708 | 0.098527 | 0.948 |
| 400  | 1.008151 | 0.004709 | 0.008151 | 0.004775 | 0.069105 | 0.946 |

# 4 Comparison with Established Models

To deepen the distribution analysis, the proposed model is compared with other established models in statistical literature. Models with similar complexity, i.e., also uniparametric models, were chosen for fair comparison of adequacy and efficiency metrics. The distributions used for validation are exponential, inverse exponential, Lindley, Rayleigh, and inverse Rayleigh.

## 4.1 Dataset
For reproducibility, the same dataset from the original article was used, which is:

|       |       |        |        |        |        |        |
|-------|-------|--------|--------|--------|--------|--------|
| 18.83 | 20.80 | 21.657 | 23.03  | 23.23  | 24.05  | 24.321 |
| 25.50 | 25.52 | 25.80  | 26.69  | 26.77  | 26.78  | 27.05  |
| 27.67 | 29.90 | 31.11  | 33.20  | 33.73  | 33.76  | 33.89  |
| 34.76 | 35.75 | 35.91  | 36.98  | 37.08  | 37.09  | 39.58  |
| 44.045| 45.29 | 45.381 | **    | **    | **    | **     |


## 4.2 Comparative Criteria

### 4.2.1 Negative Log-Likelihood
Measure of the probability that the model generated the observed data. A higher (less negative) value indicates a better fit of the model to the data.

### 4.2.2 Kolmogorov-Smirnov (K-S) Test
Maximum distance between the cumulative distribution functions of the observed data and the model. Used to test the hypothesis that the data follow a specific distribution.

### 4.2.3 Akaike Information Criterion (AIC)
Criterion based on information theory that evaluates the model's quality, balancing model fit and complexity (number of parameters). A lower AIC value is preferable.

### 4.2.4 Bayesian Information Criterion (BIC)
Similar to AIC but uses a greater penalty for the number of parameters, favoring simpler models.

### 4.2.5 Hannan-Quinn Information Criterion (HQIC)
Analogous to AIC, it follows another measure of penalization.

## 4.3 Results

### Table 5: Model Comparison - Part 1

| Distribution         | MLE -LL   | KS_stat   | KS_p           |
|----------------------|-----------|-----------|----------------|
| Exponential          | 0.032455  | 137.264447| 1.748867e-06   |
| Inverse Exponential  | 29.215334 | 137.261497| 6.155969e-07   |
| Lindley              | 0.062988  | 126.994191| 3.219054e-04   |
| Rayleigh             | 0.001000  | 118.222345| 2.651648e-03   |
| Inverse Rayleigh     | 810.503208| 118.200626| 2.015911e-03   |
| A                    | 125.662000| 107.950308| 3.543295e-01   |

### Table 6: Model Comparison - Part 2

| Distribution         | AIC       | CAIC      | BIC       | HQIC      |
|----------------------|-----------|-----------|-----------|-----------|
| Exponential          | 276.528894| 276.657926| 277.962881| 276.996338|
| Inverse Exponential  | 276.522995| 276.652027| 277.956982| 276.990439|
| Lindley              | 255.988382| 256.117414| 257.422369| 256.455826|
| Rayleigh             | 238.444691| 238.573723| 239.878678| 238.912135|
| Inverse Rayleigh     | 238.401253| 238.530285| 239.835240| 238.868697|
| A                    | 217.900616| 218.029648| 219.334603| 218.368060|

As analyzed in the base article, through classical methodology, distribution A has a lower value for all analyzed indices, suggesting better performance among the tested models.

![image](https://github.com/arthurhirez/Bayesiana/assets/109704516/a77b2a02-b915-4a3e-b63c-a71e632f1b55)

**Figure 5:** Reliability Function Plot

![image](https://github.com/arthurhirez/Bayesiana/assets/109704516/525824ac-7a6d-4e40-b627-286578da6b0a)

**Figure 6:** Hazard function Plot

![image](https://github.com/arthurhirez/Bayesiana/assets/109704516/fd93c1e5-c1b8-438f-b59c-ba466b89c8e0)

**Figure 7:** Log-Likelihood Function Plot

# 5 Study Using Bayesian Methodology

## 5.1 Markov Chain Monte Carlo (MCMC) Method

The Monte Carlo process with Markov chains involves simulating samples from posterior densities that are too complex to be obtained analytically. It starts with defining an initial point in the distribution, from which a sample is generated. In the next iteration, the result of the previous sampling process is used as the starting point, creating a Markov chain. It is assumed that Markov chains represent a "memoryless" stochastic process, depending only on the current state of the data and not on the preceding sequence of events. Convergence to a stationary and independent distribution from the initial starting point is assumed.

## 5.2 Metropolis-Hastings Algorithm

Performing stochastic simulations of a distribution without an acceptance algorithm leads to the phenomenon called "random walk," where the parameter space is explored randomly without the prospect of convergence to a desired prior. Therefore, an acceptance criterion based on the likelihood ratio is used to decide whether a new state is valid for Monte Carlo simulation.

In the implementation, at each step of MCMC, the Metropolis algorithm proposes a new value for β based on the current value. The magnitude of this proposal is influenced by the "step_scale." The likelihood ratio is then calculated, expressing which of the states is more likely to have generated the data. If greater than 1 (suggesting a better fit), it is accepted; if less than 1, it is not necessarily discarded, controlled by a random event based on the ratio. This process is repeated during the "burn-in" period, where initial values are discarded so that Markov chains "forget" the initial conditions.

## 5.3 Simulation

### 5.3.1 Uniform Prior

Due to its simple and non-informative nature, a uniform prior was initially used for the β parameter of distribution A. In this analysis, it is assumed that the choice of a uniform prior reflects a lack of prior knowledge about the problem, emphasizing the significance of the observed data in forming the posterior. The Bayesian formulation is given by:

$$ \pi(\beta|x) \propto f(x|\beta) \times \pi(\beta) \propto \frac{1}{x^2e^{\frac{1}{\beta}(1-e^{\beta x})+\beta x}} \times \frac{1}{b-a} \propto e^{\frac{1}{\beta}(1-e^{\beta x})+\beta x} $$

Values of \(a = 0\) and \(b = 100\) or \(1000\) were chosen to explore the vagueness of the prior in the study.

### 5.3.2 Gamma Prior

For comparison, a gamma distribution was also used as a prior, exploring a non-informative approach through the parameter space. The positive and unimodal nature of the distribution is highlighted. The Bayesian formulation is given by:

$$ \pi(\beta|x) \propto f(x|\beta) \times \pi(\beta) \propto \frac{1}{x^2e^{\frac{1}{\beta}(1-e^{\beta x})+\beta x}} \times \frac{b^a}{\Gamma(a)\beta^{a-1}e^{-\beta x}} \propto \beta^{a-1}e^{\frac{1}{\beta}(1-e^{\beta x})+\beta( \frac{1}{x} - b)} $$

A "flat" gamma prior was chosen, i.e., an approach that is not very informative and less restrictive, focusing on the observed data for MCMC processes. Therefore, the parameters were set to \(a = [1, 0.1]\) and \(b = [0.001, 0.0001]\). Note that the posterior in this specific case does not have an explicit kernel.

### 5.3.3 Jeffreys Prior

To obtain the Jeffreys prior - and later the estimator's variance - it is necessary to calculate the Fisher information matrix, given by the second derivative of the log-likelihood function. Due to the impossibility of finding the desired expectation analytically, numerical methods are followed. Jeffreys prior was not used in this scenario, but the model's behavior was simulated in less informative scenarios using other distributions.


### Derivatives of Log-Likelihood Function:

1. **First Derivative:**
   
$$ \frac{dl(\beta)}{d\beta} = \frac{1}{x} - \frac{1}{\beta^2} \left( \beta xe^{\beta x} - e^{\beta x} + 1 \right)^2 $$

2. **Second Derivative:**
   
$$ \frac{d^2l(\beta)}{d\beta^2} = 2e^{\beta x} \frac{\beta^3x - 2\beta^3x - 2e^{\beta x}\beta^2x^2 + e^{\beta x}\beta x^3}{3} $$

3. **Fisher Information:**
   
$$ \left[ -\frac{d^2l(\beta)}{d\beta^2} \right] = E \left[ 2e^{\beta x} \frac{\beta^3x - 2\beta^3x - 2e^{\beta x}\beta^2x^2 + e^{\beta x}\beta x^3}{3} \right] $$

   Due to the impossibility of finding the desired expectation analytically, numerical methods are used.

Considering this, the Jeffreys prior was not used. However, the behavior of the model was simulated in less informative scenarios using other distributions.

## 5.4 Results and Discussions

Experiments were conducted for $(\beta = [0.125, 0.6, 1]\)$ to compare the Bayesian approach to the classical one. The goal, once again, is to calculate the same bias and variance metrics, this time obtained through Monte Carlo simulations.

**Table 7:** Simulation Results with Gamma Prior (a = 1, b = 0.001) for β = 0.125
| Size | Mean     | Var      | Bias     | MSE      | RMSE     | PC    |
|------|----------|----------|----------|----------|----------|-------|
| 25   | 0.255756 | 0.028168 | 0.130756 | 0.045265 | 0.212756 | 0.929 |
| 50   | 0.195230 | 0.013095 | 0.070230 | 0.018027 | 0.134266 | 0.937 |
| 100  | 0.153991 | 0.005144 | 0.028991 | 0.005985 | 0.077362 | 0.953 |
| 200  | 0.136739 | 0.002824 | 0.011739 | 0.002961 | 0.054419 | 0.936 |
| 400  | 0.128950 | 0.001574 | 0.003950 | 0.001590 | 0.039873 | 0.927 |

**Table 8:** Simulation Results with Uniform Prior (a = 0, b = 1000.0) for β = 0.125
| Size | Mean     | Var      | Bias     | MSE      | RMSE     | PC    |
|------|----------|----------|----------|----------|----------|-------|
| 25   | 0.255122 | 0.028153 | 0.130122 | 0.045085 | 0.212332 | 0.931 |
| 50   | 0.194912 | 0.013087 | 0.069912 | 0.017974 | 0.134069 | 0.939 |
| 100  | 0.153862 | 0.005150 | 0.028862 | 0.005983 | 0.077350 | 0.949 |
| 200  | 0.136699 | 0.002827 | 0.011699 | 0.002964 | 0.054441 | 0.934 |
| 400  | 0.128953 | 0.001574 | 0.003953 | 0.001590 | 0.039874 | 0.928 |

**Table 9:** Simulation Results with Gamma Prior (a = 1, b = 0.001) for β = 0.6
| Size | Mean     | Var      | Bias     | MSE      | RMSE     | PC    |
|------|----------|----------|----------|----------|----------|-------|
| 25   | 0.666222 | 0.066290 | 0.066222 | 0.070676 | 0.265849 | 0.929 |
| 50   | 0.629221 | 0.033831 | 0.029221 | 0.034685 | 0.186240 | 0.934 |
| 100  | 0.607921 | 0.013651 | 0.007921 | 0.013713 | 0.117104 | 0.942 |
| 200  | 0.604802 | 0.006706 | 0.004802 | 0.006729 | 0.082032 | 0.949 |
| 400  | 0.602869 | 0.003295 | 0.002869 | 0.003303 | 0.057473 | 0.941 |

**Table 10:** Simulation Results with Uniform Prior (a = 0, b = 1000.0) for β = 0.6
| Size | Mean     | Var      | Bias     | MSE      | RMSE     | PC    |
|------|----------|----------|----------|----------|----------|-------|
| 25   | 0.666282 | 0.066443 | 0.066282 | 0.070837 | 0.266151 | 0.929 |
| 50   | 0.629130 | 0.033829 | 0.029130 | 0.034678 | 0.186219 | 0.933 |
| 100  | 0.607899 | 0.013640 | 0.007899 | 0.013702 | 0.117055 | 0.944 |
| 200  | 0.604710 | 0.006701 | 0.004710 | 0.006723 | 0.081992 | 0.948 |
| 400  | 0.602876 | 0.003297 | 0.002876 | 0.003305 | 0.057492 | 0.942 |

**Table 11:** Simulation Results with Gamma Prior (a = 1, b = 0.001) for β = 1.0
| Size | Mean     | Var      | Bias     | MSE      | RMSE     | PC    |
|------|----------|----------|----------|----------|----------|-------|
| 25   | 1.056040 | 0.094632 | 0.056040 | 0.097773 | 0.312686 | 0.940 |
| 50   | 1.028186 | 0.047459 | 0.028186 | 0.048254 | 0.219667 | 0.931 |
| 100  | 1.008461 | 0.019106 | 0.008461 | 0.019178 | 0.138483 | 0.960 |
| 200  | 1.005019 | 0.009451 | 0.005019 | 0.009477 | 0.097348 | 0.947 |
| 400  | 1.002835 | 0.004702 | 0.002835 | 0.004710 | 0.068627 | 0.945 |

**Table 12:** Simulation Results with Uniform Prior (a = 0, b = 1000.0) for β = 1.0
| Size | Mean     | Var      | Bias     | MSE      | RMSE     | PC    |
|------|----------|----------|----------|----------|----------|-------|
| 25   | 1.055836 | 0.094798 | 0.055836 | 0.097915 | 0.312914 | 0.939 |
| 50   | 1.028129 | 0.047484 | 0.028129 | 0.048275 | 0.219717 | 0.930 |
| 100  | 1.008404 | 0.019111 | 0.008404 | 0.019182 | 0.138498 | 0.955 |
| 200  | 1.005010 | 0.009459 | 0.005010 | 0.009484 | 0.097387 | 0.947 |
| 400  | 1.002869 | 0.004708 | 0.002869 | 0.004716 | 0.068673 | 0.948 |

## 5.5 Bayesian Approach to the Original Data

It is observed that the coverage probability is calculated differently in the Bayesian paradigm. In contrast to the classical methodology, where the coverage probability represents the frequency with which the true parameter value belongs to the confidence interval calculated given a large number of repetitions, here PC is the probability that the credibility interval contains the true parameter value, given the posterior distribution.

### 5.5 Data from the Article

With the intention of repeating the steps taken in the original article, the same dataset is analyzed from a Bayesian perspective for the presented priors:

#### 5.5.1 Prior Gamma

The histogram below displays the distribution of samples, providing an insight into the posterior distribution of the parameter. The central tendency of the samples around 125 is noticeable, indicating good results.

![image](https://github.com/arthurhirez/Bayesiana/assets/109704516/3ea260d0-730c-42a6-8caa-cae643d5d8db)

**Figure 8:** Plot of the estimated beta distribution for the Gamma prior

The trace plot shows the trajectories of the samples over the iterations of the MCMC for each chain. It allows visualizing the convergence of the simulation over time. In this example, adequate convergence is indicated by chains that seem to mix well, without discernible trends or patterns.

![image](https://github.com/arthurhirez/Bayesiana/assets/109704516/4e5779f3-6589-4455-b0b5-2cc33a758b68)

**Figure 9:** Plot of the beta value over the Markov chain for the Gamma prior

In the autocorrelation plot, a tendency to reduce autocorrelation over iterations is observed, revealing stabilization close to zero. This indicates that after the original burn-in period, the samples become independent.

![image](https://github.com/arthurhirez/Bayesiana/assets/109704516/a0ac6cf2-2bf0-44b9-83f9-b63e66ed54be)

**Figure 10:** Autocorrelation plot

##### Summary Statistics - Part 1

| Parameter | Mean   | SD    | HDI 3%  | HDI 97%  |
|-----------|--------|-------|---------|----------|
| beta      | 125.208| 5.031 | 115.859 | 134.685  |

##### Summary Statistics - Part 2

| Parameter | MCSE Mean | MCSE SD | ESS Bulk |  ESS Tail | R-hat |
|-----------|-----------|---------|----------|-----------|-------|
| beta      | 0.078     | 0.055   | 4224.0   | 4553.0    | 1.0   |

#### 5.5.2 Uniform Prior

Analogous to the flat Gamma prior, the uniform prior shows its central tendency around 125, suggesting adequacy.

#### Summary Statistics - Part 1

| Parameter | Mean   | SD    | HDI 3%  | HDI 97% |
|-----------|--------|-------|---------|---------|
| beta      | 125.147| 5.035 | 115.455 | 134.27  |

#### Summary Statistics - Part 2

| Parameter | MCSE Mean | MCSE SD | ESS Bulk |  ESS Tail | R-hat |
|-----------|-----------|---------|----------|-----------|-------|
| beta      | 0.083     | 0.059   | 3683.0   | 3981.0    | 1.0   |

![image](https://github.com/arthurhirez/Bayesiana/assets/109704516/a73a5908-7d68-4437-964c-9084c7864e26)

**Figure 11:** Plot of the estimated beta distribution for the uniform prior

The trajectory of the samples also does not reveal patterns, cycles, or trends.

![image](https://github.com/arthurhirez/Bayesiana/assets/109704516/3b7a1c49-6b37-40de-a06d-c715281e9c48)

**Figure 12:** Plot of the beta value over the Markov chain for the uniform prior

![image](https://github.com/arthurhirez/Bayesiana/assets/109704516/48f65ec5-2c55-4e6e-96b4-35466370cd2b)

**Figure 13:** Autocorrelation plot



## 6. Conclusion

Initially, it is important to note that the selection of the prior and the parameters for the prior is a crucial step, as this choice directly impacts the behavior of the Bayesian model. Therefore, the use of uninformative priors characterized by high variance is opted for due to the absence of an expert in the field.

### 6.1 Simulations

When analyzing the simulations, it becomes evident that the lower the value of the beta parameter, the worse the convergence of the Bayesian model concerning the log-likelihood calculation, as illustrated in tables 1, 7, and 8. On the other hand, for higher values of beta, the parameters estimated by the model converge more quickly to the true parameter value, as presented in tables 3, 9, and 10. However, it is worth noting that changes in the values of the prior parameters influence these results, indicating that more accurate results can be obtained with smaller betas.

### 6.2 Data from the Article

By using the same priors to estimate the beta value, we obtained the following results: Prior Gamma: 125.208, Prior Uniform: 125.147. Both values are close; however, they are lower than the predicted value using Maximum Likelihood Estimation (MLE). This observation reinforces, as highlighted in the Bayesian simulation tables, that both uninformative priors exhibit similar behaviors. It is notable that both priors converged rapidly, as evidenced in figures 9 and 12, and generated several poorly autocorrelated estimators, as demonstrated in figures 10 and 13.

### 6.3 Final Considerations

As observed, Bayesian models provide robust and meaningful results. However, it is crucial to make the appropriate choice of priors and parameters for each problem. Additionally, it is important to note that these models tend to be computationally more expensive compared to classical methods. Therefore, mastering both tools and using them at the most appropriate moments becomes extremely important for an effective and efficient approach.

