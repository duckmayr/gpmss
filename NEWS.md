# gpmss 0.1.1

## Bug fixes

* Gracefully handle factor levels not present in training data for `predict()`.
  When `GPModel` class objects are created, factor variables are one-hot
  encoded, using all factor levels, meaning some factor levels observed in
  testing data may not be observed in training data. This possibility is now
  handled without error.

# gpmss 0.1.0

* Initial release
* Regression and binary classification are supported
* Mean zero and a linear mean function are available
* The isometric squared exponential and the automatic relevance determination
  squared exponential kernels are available for covariance functions
