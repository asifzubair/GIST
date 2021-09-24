//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data
data
{
  int<lower=0> numGenes; // The number of "genes" (1000)
  int<lower=0> numCellTypes; // number of cell types (2)
  vector[numGenes] exprMixVec; // the data

  // the matrix of signature (numGenes, numCellTypes)
  matrix[numGenes, numCellTypes] sigMat;
}

transformed data
{
  vector<lower=0>[numCellTypes] alpha;
  for (i in 1:numCellTypes)
  {
    // This is an uniform prior over the dirichlet (beta distribution if numCellTypes = 2).
    // Vector of 1s is uniform in the case of  dirichlet distrubution.
    alpha[i] = 1;
  }
}

// The parameters accepted by the model.
parameters
{
  // Syntax is "simplex[dimensions_of_simplexes] vectorOfSimplexesName[length_of_vector]".
  simplex[numCellTypes] estimatedProportionsVecSimp;

  real<lower=1> nu;

  // Variance parameters, estimated below.
  real<lower=0> sigma;
  real beta0;
}

// The model to be estimated.
model
{
  // alpha is a vector 1's, meaning this dirichlet represents a uniform prior.
  estimatedProportionsVecSimp ~ dirichlet(alpha);

  // Note that the gamma distribution is in terms of shape and rate
  // i.e. y ~ gamma(alpha, beta)
  // this means that for our case: E(nu) = 2/0.01 = 200 !
  nu ~ gamma(2, 0.1);

  exprMixVec ~ student_t(nu, beta0 + sigMat * estimatedProportionsVecSimp, sigma);

}
