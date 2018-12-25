/*****/
// Normal distribution
/*****/

var normalMean = new unboundedParameter(
	"mu",
	"\\mu",
	function() { return [-5,5]; },
	false,
	0,
	"continuous",
	"normal",
	"location"
);

var normalVariance = new positiveParameter(
	"sigma2",
	"\\sigma^2",
	function() { return [0.1,10]; },
	true,
	1,
	"continuous",
	"inversegamma",
	"variance"
);

var normalStandardDeviation = new positiveParameter(
	"sigma",
	"\\sigma",
	function() { return [0.1,10]; },
	true,
	1,
	"continuous",
	null,
	"standard deviation"
);


var normalPrecision = new positiveParameter(
	"tau",
	"\\tau",
	function() { return [0.1,10]; },
	true,
	1,
	"continuous",
	"gamma",
	"precision"
);

var normalMeanVariance = new distributionParametrization(
	"mean/variance",
	[normalMean, normalVariance],
	function(mu,sig2) { return [Number.NEGATIVE_INFINITY,Number.POSITIVE_INFINITY]; },
	function(mu,sig2) { return [-5,5]; },
	"(-\\infty,\\infty)",
	function(x, mu, sig2) { return jStat.normal.pdf(x, mu,Math.sqrt(sig2)); },
	"\\left(2\\pi\\sigma^2\\right)^{-\\frac{1}{2}}\\exp\\left\\{-\\frac{1}{2\\sigma^2}\\left(x-\\mu\\right)^2\\right\\}",
	function(x, mu, sig2) { return jStat.normal.cdf(x, mu,Math.sqrt(sig2)); },
	"normalinversegamma",
	function(mu,sig2){ return sig2>0;},
	{
		mean: {
			fun: function(mu, sig2) { return jStat.normal.mean(mu, Math.sqrt(sig2)); },
			display: "\\mu"
			},
		variance: {
			fun: function(mu, sig2) { return jStat.normal.variance(mu, Math.sqrt(sig2)); },
			display: "\\sigma^2"
			},
		median: {
			fun: function(mu, sig2) { return jStat.normal.median(mu, Math.sqrt(sig2)); },
			display: "\\mu"
			},
		mode: {
			fun: function(mu, sig2) { return jStat.normal.mode(mu, Math.sqrt(sig2)); },
			display: "\\mu"
			}

	},
	null
);

var normalMeanPrecision = new distributionParametrization(
	"mean/precision",
	[normalMean, normalPrecision],
	function(mu,tau) { return [Number.NEGATIVE_INFINITY,Number.POSITIVE_INFINITY]; },
	function(mu,tau) { return [-5,5]; },
	"(-\\infty,\\infty)",
	function(x,mu,tau){ return jStat.normal.pdf(x,mu,Math.sqrt(1/tau)); },
	"\\left(\\frac{\\tau}{2\\pi}\\right)^{\\frac{1}{2}}\\exp\\left\\{-\\frac{\\tau}{2}\\left(x-\\mu\\right)^2\\right\\}",
	function(x,mu,tau){ return jStat.normal.cdf(x,mu,Math.sqrt(1/tau)); },
	"normalgamma",
	function(mu,tau){ return tau>0;},
	{
		mean: {
			fun: function(mu, tau) { return jStat.normal.mean(mu, 1/Math.sqrt(tau)); },
			display: "\\mu"
			},
		variance: {
			fun: function(mu, tau) { return jStat.normal.variance(mu, 1/Math.sqrt(tau)); },
			display: "\\frac{1}{\\tau}"
			},
		median: {
			fun: function(mu, tau) { return jStat.normal.median(mu, 1/Math.sqrt(tau)); },
			display: "\\mu"
			},
		mode: {
			fun: function(mu, tau) { return jStat.normal.mode(mu, 1/Math.sqrt(tau)); },
			display: "\\mu"
			}
	},
	null

);

var normalMeanStandardDeviation = new distributionParametrization(
	"mean/standard deviation",
	[normalMean, normalStandardDeviation],
	function(mu,sig2) { return [Number.NEGATIVE_INFINITY,Number.POSITIVE_INFINITY]; },
	function(mu,sig2) { return [-5,5]; },
	"(-\\infty,\\infty)",
	function(x, mu, sig) { return jStat.normal.pdf(x, mu, sig); },
	//1/sqrt(2*)
	"\\left(2\\pi\\sigma^2\\right)^{-\\frac{1}{2}}\\exp\\left\\{-\\frac{1}{2\\sigma^2}\\left(x-\\mu\\right)^2\\right\\}",
	function(x, mu, sig) { return jStat.normal.cdf(x, mu, sig); },
	null,
	function(mu,sig){ return sig>0;},
	{
		mean: {
			fun: function(mu, sig) { return jStat.normal.mean(mu, sig); },
			display: "\\mu"
			},
		variance: {
			fun: function(mu, sig) { return jStat.normal.variance(mu, sig); },
			display: "\\sigma"
			},
		median: {
			fun: function(mu, sig) { return jStat.normal.median(mu, sig); },
			display: "\\mu"
			},
		mode: {
			fun: function(mu, sig) { return jStat.normal.mode(mu, sig); },
			display: "\\mu"
			}

	},
	null
);


distributions["normal"] = new distribution(
	"normal",
	"Normal/Gaussian",
	"continuous",
	[ normalMeanVariance, normalMeanPrecision, normalMeanStandardDeviation],
	null,
	{}
);






/*****/
// Gamma distribution
/*****/

var gammaShape = new positiveParameter(
	"k",
	"k",
	function() { return [.1,10]; },
	true,
	1,
	"continuous",
	null,
	"shape"
);

var gammaScale = new positiveParameter(
	"s",
	"s",
	function() { return [0.1,10]; },
	true,
	1,
	"continuous",
	"inversegamma",
	"scale"
);

var gammaRate = new positiveParameter(
	"theta",
	"\\theta",
	function() { return [0.1,10]; },
	true,
	1,
	"continuous",
	"gamma",
	"rate"
);

var gammaShapeScale = new distributionParametrization(
	"shape/scale",
	[gammaShape, gammaScale],
	function(k,s) { return [0,Number.POSITIVE_INFINITY]; },
	function(k,s) { return [.01,10]; },
	"(0,\\infty)",
	jStat.gamma.pdf,
	"\\frac{1}{\\Gamma(k)s^{k}} x^{k - 1} \\exp\\left\\{-\\frac{x}{s}\\right\\}",
	jStat.gamma.cdf,
	null,
	function(k,s){ return s>0 && k>0;},
	{
		mean: {
			fun: jStat.gamma.mean,
			display: "ks"
			},
		variance: {
			fun: jStat.gamma.variance,
			display: "ks^2"
			},
		mode: {
			fun: jStat.gamma.mode,
			display: "(k-1)s, k\\geq 1"
			}
	},
	null

);

var gammaShapeRate = new distributionParametrization(
	"shape/rate",
	[gammaShape, gammaRate],
	function(k,theta) { return [0,Number.POSITIVE_INFINITY]; },
	function(k,theta) { return [.01,10]; },
	"(0,\\infty)",
	function(x, k, theta) { return jStat.gamma.pdf(x, k, 1/theta); },
	"\\frac{\\theta^k}{\\Gamma(k)} x^{k - 1} \\exp\\left\\{-\\theta x\\right\\}",
	function(x, k, theta) { return jStat.gamma.cdf(x, k, 1/theta); },
	null,
	function(k,theta){ return theta>0 && k>0;},
	{
		mean: {
			fun: function(k,theta) { return jStat.gamma.mean(k, 1/theta); },
			display: "\\frac{k}{\\theta}"
			},
		variance: {
			fun: function(k,theta) { return jStat.gamma.variance(k, 1/theta); },
			display: "\\frac{k}{\\theta^2}"
			},
		mode: {
			fun: function(k,theta) { return jStat.gamma.mode(k, 1/theta); },
			display: "\\frac{k-1}{\\theta}, k\\geq 1"
			}
	},
	null

);
distributions["gamma"] = new distribution(
	"gamma",
	"Gamma",
	"continuous",
	[ gammaShapeScale, gammaShapeRate ],
	null,
	{
		// link:"http://en.wikipedia.org/wiki/Gamma_distribution"
	}
);




/*****/
// Binomial distribution
/*****/

var binomialN = new positiveParameter(
	"N",
	"N",
	function() { return [1,100]; },
	false,
	20,
	"discrete",
	null,
	"sample size"
);

var binomialp = new distributionParameter(
	"p",
	"p",
	function() { return [0,1]; },
	function() { return [0,1]; },
	false,
	.5,
	"continuous",
	"beta",
	"probability of success"
);

var binomialNp = new distributionParametrization(
	"probability",
	[binomialN, binomialp],
	function(N, p) { return [0,N]; },
	function(N, p) { return [0,N]; },
	"[0,N]",
	jStat.binomial.pdf,
	"\\binom{N}{x} p^x (1-p)^{N-x}",
	jStat.binomial.cdf,
	"beta",
	function(N,p){ return N>0 && is_int(N) && p<1 && p>0;},
	{
		mean: { 
			fun: jStat.binomial.mean, 
			display: "Np" 
			},
		variance: {
			fun: jStat.binomial.variance, 
			display: "Np(1-p)" 
			}
	},
	null

);

distributions["binomial"] = new distribution(
	"binomial",
	"Binomial",
	"discrete",
	[ binomialNp ],
	null,
	{

	}

);


/*****/
// Poisson distribution
/*****/

var poissonRate = new positiveParameter(
	"lambda",
	"\\lambda",
	function() { return [1/10,10]; },
	true,
	1,
	"continuous",
	null,
	"rate"
);

var poissonRate = new distributionParametrization(
	"rate",
	[ poissonRate ],
	function(lambda) { return [0,Number.POSITIVE_INFINITY]; },
	function(lambda) { return [0,25]; },
	"[0,\\infty)",
	jStat.poisson.pdf,
	"\\frac{\\lambda^x}{x!} \\exp\\left\\{-\\lambda\\right\\}",
	jStat.poisson.cdf,
	"gamma",
	function(lambda){ return lambda>0;},
	{
		mean: {
			fun: jStat.poisson.mean,
			display: "\\lambda"
			},
		variance: {
			fun: jStat.poisson.variance,
			display: "\\lambda"
			}
	},
	null
);

distributions["poisson"] = new distribution(
	"poisson",
	"Poisson",
	"discrete",
	[ poissonRate ],
	null,
	{

	}
);



/*****/
// Negative Binomial distribution
/*****/

var negbinomialr = new positiveParameter(
	"r",
	"r",
	function() { return [1,50]; },
	false,
	20,
	"discrete",
	null,
	"number of successes required"
);


var negbinomialp = new distributionParameter(
	"p",
	"p",
	function() { return [0,1]; },
	function() { return [0,1]; },
	false,
	.2,
	"continuous",
	"beta",
	"probability of failure"
);
var binomialNp = new distributionParametrization(
    "probability",
    [binomialN, binomialp],
    function(N, p) { return [0,N]; },
    function(N, p) { return [0,N]; },
    "[0,N]",
    jStat.binomial.pdf,
    "\\binom{N}{x} p^x (1-p)^{N-x}",
    jStat.binomial.cdf,
    "beta",
    function(N,p){ return N>0 && is_int(N) && p<1 && p>0;},
    {
        mean: {
            fun: jStat.binomial.mean,
            display: "Np"
        },
        variance: {
            fun: jStat.binomial.variance,
            display: "Np(1-p)"
        }
    },
    null

);

var negbinomialrp = new distributionParametrization(
	"probability",
	[negbinomialr, negbinomialp],
	function(r, p) { return [0,Number.POSITIVE_INFINITY]; },
	function(r, p) { return [0,50]; },
	"[0,\\infty)",
	jStat.negbin.pdf,
	"\\binom{x+r-1}{x} p^x (1-p)^{r}",
	jStat.negbin.cdf,
	"beta",
	function(r,p){ return r>0 && is_int(r) && p<1 && p>0;},
	{
		mean: {
			fun: function(r,p){ return p*r/(1-p);},
			display: "\\frac{rp}{1-p}"
			},
		variance: {
			fun: function(r,p){ return p*r/(1-p)^2;},
			display: "\\frac{rp}{(1-p)^2}"
			}
	},
	null

);

distributions["negativebinomial"] = new distribution(
	"negativebinomial",
	"Negative binomial",
	"discrete",
	[ negbinomialrp ],
	null,
	{

	}
);
// mean: {
//     fun: function(r,p){ return p*r/(1-p);},
//     display: "\\frac{rp}{1-p}"
// },
// variance: {
//     fun: function(r,p){ return p*r/(1-p)^2;},
//     display: "\\frac{rp}{(1-p)^2}"
// }
