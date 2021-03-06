## -----------------------------------------------------------------------------
    library(lab4)

## -----------------------------------------------------------------------------
linreg_mod <- linreg$new(Form = Petal.Length~Species, data=iris)

## -----------------------------------------------------------------------------
linreg_mod$print()

## -----------------------------------------------------------------------------
linreg_mod$plot()

## -----------------------------------------------------------------------------
head(linreg_mod$resid())

## -----------------------------------------------------------------------------
head(linreg_mod$pred())

## -----------------------------------------------------------------------------
linreg_mod$coef()

## -----------------------------------------------------------------------------
linreg_mod$summary()

