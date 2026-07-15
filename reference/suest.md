# Combine Two Fitted Models with Seemingly Unrelated Estimation

Combines two separately fitted models and constructs a joint
model-robust covariance matrix from their observation-level score
contributions.

## Usage

``` r
suest(model1, model2, model_names = NULL)

# S3 method for class 'suest_model'
coef(object, ...)

# S3 method for class 'suest_model'
vcov(object, ...)

# S3 method for class 'suest_model'
nobs(object, ...)

# S3 method for class 'suest_model'
print(x, ...)
```

## Arguments

- model1, model2:

  Two supported fitted model objects.

- model_names:

  Optional character vector containing two display names.

- object, x:

  A `suest_model` object.

- ...:

  Additional arguments.

## Value

An object of class `"suest_model"` containing the component models, a
joint coefficient vector, a joint model-robust covariance matrix, and
sample-alignment information.

## Details

Exactly two models are supported. Models may use identical, partially
overlapping, or disjoint samples. Supported model classes include
linear, binary logit and probit, Poisson, negative binomial, ordered
logit and probit, and multinomial logit models.

The returned object works with the marginaleffects package.

## Examples

``` r
dat <- mtcars
dat$am <- factor(dat$am)
model1 <- glm(am ~ wt, family = binomial(), data = dat)
model2 <- glm(am ~ wt + hp, family = binomial(), data = dat)
fit <- suest(model1, model2, model_names = c("Base", "Adjusted"))
fit
#> Seemingly Unrelated Estimation
#> Models: Base + Adjusted 
#> Model types: Base=logit, Adjusted=logit 
#> Comparison scale: predicted probabilities 
#> Observations: Base=32, Adjusted=32 
#> Overlapping observations: 32 
#> Union observations: 32 
#> Parameters: 5 
```
