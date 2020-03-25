import numpy as np
import scipy.optimize

def softmax(X):
    '''
    Arg: numpy.ndarray. shape=(n_samples, n_classes)
      Score.
      If the input is a vector, it is assumed to be the score of the "positive" class.
    
    Returns: numpy.ndarray. shape=(n_samples, n_classes)
      Output shape is always the same as the input shape.
      If the input is a vector, returns the probability of the "positive" class.
    '''
    if X.ndim == 1 or (X.ndim == 2 and X.shape[1] == 1):
        return 1 / (1 + np.exp(-X))
    exps = np.exp(X - np.max(X, axis=1, keepdims=True))
    return exps / np.sum(exps, axis=1, keepdims=True)

def cross_entropy(pred, y, validate=False, sample_weight=None):
    '''
    Cross entropy loss: H(pred, y) = - 1/n_samples * sum_{i in samples} sum_{j in classes} y[i,j] * log(pred[i,j])
    
    Args
    - pred: numpy.ndarray. shape=(n_samples, n_classes)
        Distribution over classes
    - y: numpy.ndarray. shape=(n_samples, n_classes)
        Target probability distribution over classes
    - sample_weight: numpy.ndarray. shape=(n_samples,)
        Sample weights
    
    Returns: float
    '''
    if y.ndim == 1:
        y = np.vstack((1-y, y)).T
    elif y.ndim == 2 and y.shape[1] == 1:
        y = np.hstack((1-y, y))
    
    if pred.ndim == 1:
        pred = np.vstack((1-pred, pred)).T
    elif pred.ndim == 2 and pred.shape[1] == 1:
        pred = np.hstack((1-pred, pred))

    assert pred.shape == y.shape

    if validate:
        assert np.all(np.sum(y, axis=1) == 1), '`y` probabilities must sum to 1 over all classes.'
        assert np.all(np.sum(pred, axis=1) == 1), '`pred` probabilities must sum to 1 over all classes.'
    n_samples = y.shape[0]
    with np.errstate(divide='ignore'):
        log_pred = np.log(pred)
    log_pred[~np.isfinite(log_pred)] = 0
    if sample_weight is not None:
        assert len(sample_weight) == pred.shape[0]
        if y.ndim == 2 and sample_weight.ndim == 1:
            sample_weight = np.expand_dims(sample_weight, axis=1)
        loss = - np.sum(y * log_pred * sample_weight) / n_samples
    else:
        loss = - np.sum(y * log_pred) / n_samples
    return loss

class logistic_regression:
    '''
    Logistic regression with soft labels.
    Largely modeled off of scikit-learn's LogisticRegression implementation.
    '''
    def __init__(self, alpha=0.001, beta=0, bias=True, method='L-BFGS-B', optim_kwargs=None, rng=None):
        if rng is None:
            rng = np.random.default_rng()
        if optim_kwargs is None:
            optim_kwargs = {}
        self.alpha = alpha
        self.beta = beta
        self.bias = bias
        self.method = method
        self.optim_kwargs = optim_kwargs
        self.rng = rng

    def fit(self, X, y, sample_weight=None, weights0=None, **kwargs):
        assert X.shape[0] == y.shape[0]

        if self.bias:
            X = np.hstack((np.ones((X.shape[0], 1)), X))
        n_features = X.shape[1]
        if y.ndim == 1 or (y.ndim == 2 and y.shape[1] == 1):
            n_classes = 2
        else:
            n_classes = y.shape[1]
        
        if weights0 is None:
            weights0 = self.rng.random((n_classes - 1, n_features))

        self.result_ = scipy.optimize.minimize(
            self._loss,
            weights0,
            args=(X, y, sample_weight),
            method=self.method,
            jac=True,
            options=self.optim_kwargs,
            **kwargs)

        self.w_ = self.result_.x.reshape(n_classes - 1, n_features)
        if self.bias:
            self.intercept_ = self.w_[:, 0]
            self.coef_ = self.w_[:, 1:]
        else:
            self.coef_ = self.w_.copy()
            self.intercept_ = 0

    def _loss(self, weights, X, y, sample_weight):
        '''
        Args
        - weights: numpy.ndarray. shape=(n_classes * n_features,)
        - X: numpy.ndarray. shape=(n_samples, n_features)
        - y: numpy.ndarray. shape=(n_samples, n_classes)
        - sample_weight: numpy.ndarray. shape=(n_samples,)
        '''
        if sample_weight is None:
            sample_weight = 1
        elif y.ndim == 2 and sample_weight.ndim == 1:
            sample_weight = np.expand_dims(sample_weight, axis=1)
        if y.ndim == 1:
            y = np.vstack((1-y, y)).T

        n_features = X.shape[1]

        # the factor of 0.5 for the L2 term follows scikit-learn's implementation
        loss_reg = self.alpha * 0.5 * np.linalg.norm(weights) + \
                   self.beta * np.linalg.norm(weights, ord=1)

        weights = np.vstack((np.zeros((1, n_features)), weights.reshape(-1, n_features)))
        # softmax is over-parameterized; let first class always have score 0
        score = X @ weights.T
        
        logp = score - scipy.special.logsumexp(score, axis=1, keepdims=True)
        loss_ce = -np.sum(sample_weight * y * logp)
        # alternatively
        # p = softmax(score)
        # loss_ce = cross_entropy(pred, y, sample_weight)
        
        diff = sample_weight * (np.exp(logp) - y)
        grad = diff.T @ X
        grad += self.alpha * weights
        grad += self.beta * np.sign(weights)
        grad = grad[1, :] - grad[0, :]

        return loss_ce + loss_reg, grad

    def predict(self, X, squeeze=True):
        if self.bias:
            X = np.hstack((np.ones((X.shape[0], 1)), X))
        n_samples = X.shape[0]
        pred = softmax(np.hstack((np.zeros((n_samples, 1)), X @ self.w_.T)))
        if squeeze:
            return np.squeeze(pred[:, 1:])
        return pred

    def get_params(self, deep=True):
        return {'alpha': self.alpha, 'beta': self.beta}

    def set_params(self, **parameters):
        for parameter, value in parameters.items():
            setattr(self, parameter, value)
        return self