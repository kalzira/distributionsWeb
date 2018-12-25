/**
 * jStat */
this.j$ = this.jStat = (function( Math, undefined ) {

	// for quick reference
var slice = Array.prototype.slice,
	toString = Object.prototype.toString,

	// ascending/descending functions for sort
	ascNum = function( a, b ) { return a - b; },

	// calculate correction for IEEE error
	calcRdx = function( n, m ) {
		var val = n > m ? n : m;
		return Math.pow( 10, 17 - ~~( Math.log((( val > 0 ) ? val : -val )) * Math.LOG10E ));
	},

	// test if array
	isArray = Array.isArray || function( arg ) {
		return toString.call( arg ) === "[object Array]";
	},

	// test if function
	isFunction = function( arg ) {
		return toString.call( arg ) === "[object Function]";
	};

// global function
function jStat() {
	return new jStat.fn.init( arguments );
}

// extend jStat prototype
jStat.fn = jStat.prototype = {
	constructor : jStat,
	init : function( args ) {
		// if first argument is an array, must be vector or matrix
		if ( isArray( args[0] )) {
			// check if matrix
			if ( isArray( args[0][0] )) {
				for ( var i = 0; i < args[0].length; i++ ) {
					this[i] = args[0][i];
				}
				this.length = args[0].length;
			// so must be vector
			} else {
				this[0] = args[0];
				this.length = 1;
			}
		// if first argument is number, assume creation of sequence
		} else if ( !isNaN( args[0] )) {
			this[0] = jStat.seq.apply( null, args );
			this.length = 1;
		}
		return this;
	},

	// default length
	length : 0,

	// return clean array
	toArray : function() {
		return ( this.length > 1 ) ?
			slice.call( this )
		: slice.call( this )[0];
	},

	// only to be used internally
	push : [].push,
	sort : [].sort,
	splice : [].splice
};

// for later instantiation
jStat.fn.init.prototype = jStat.fn;

// utility functions
jStat.utils = {
	calcRdx : calcRdx,
	isArray : isArray,
	isFunction : isFunction
};

// create method for easy extension
jStat.extend = function( obj ) {
	var args = slice.call( arguments ),
		i = 1, j;
	if ( args.length === 1 ) {
		for ( j in obj ) {
			jStat[j] = obj[j];
		}
		return this;
	}
	for ( ; i < args.length; i++ ) {
		for ( j in args[i] ) obj[j] = args[i][j];
	}
	return obj;
};

// static methods
jStat.extend({

	// transpose a matrix or array
	transpose : function( arr ) {
		if ( !isArray( arr[0] )) arr = [ arr ];
		var rows = arr.length,
			cols = arr[0].length,
			obj = [],
			i = 0, j;
		for ( ; i < cols; i++ ) {
			obj.push([]);
			for ( j = 0; j < rows; j++ ) {
				obj[i].push( arr[j][i] );
			}
		}
		return obj;
	},


	identity : function( rows, cols ) {
		if ( isNaN( cols )) cols = rows;
		return jStat.create( rows, cols, function( i, j ) { return ( i === j ) ? 1 : 0; });
	},

	// generate sequence
	seq : function( min, max, length, func ) {
		if ( !isFunction( func )) func = false;
		var arr = [],
			hival = calcRdx( min, max ),
			step = ( max * hival - min * hival ) / (( length - 1 ) * hival ),
			current = min,
			cnt = 0;
		// current is assigned using a technique to compensate for IEEE error
		for ( ; current <= max; cnt++, current = ( min * hival + step * hival * cnt ) / hival )
			arr.push(( func ? func( current, cnt ) : current ));
		return arr;
	},

	// add a vector/matrix to a vector/matrix or scalar
	add : function( arr, arg ) {
		// check if arg is a vector or scalar
		if ( isArray( arg )) {
			if ( !isArray( arg[0] )) arg = [ arg ];
			return jStat.map( arr, function( value, row, col ) { return value + arg[row][col]; });
		}
		return jStat.map( arr, function ( value ) { return value + arg; });
	},




	// mean value of an array
	mean : function( arr ) {
		return jStat.sum( arr ) / arr.length;
	},

	// median of an array
	median : function( arr ) {
		var arrlen = arr.length,
			_arr = arr.slice().sort( ascNum );
		// check if array is even or odd, then return the appropriate
		return !( arrlen & 1 )
			? ( _arr[( arrlen / 2 ) - 1 ] + _arr[( arrlen / 2 )]) / 2
		: _arr[( arrlen / 2 ) | 0 ];
	},
	mode : function( arr ) {
		var arrLen = arr.length,
			_arr = arr.slice().sort( ascNum ),
			count = 1,
			maxCount = 0,
			numMaxCount = 0,
			i = 0,
			maxNum;
		for ( ; i < arrLen; i++ ) {
			if ( _arr[ i ] === _arr[ i + 1 ] ) {
				count++;
			} else {
				if ( count > maxCount ) {
					maxNum = _arr[i];
					maxCount = count;
					count = 1;
					numMaxCount = 0;
				} else {
					// are there multiple max counts
					if ( count === maxCount ) {
						numMaxCount++;
					// count is less than max count, so reset values
					} else {
						count = 1;
					}
				}
			}
		}
		return ( numMaxCount === 0 ) ? maxNum : false;
	},

	// range of an array
	range : function( arr ) {
		var _arr = arr.slice().sort( ascNum );
		return _arr[ _arr.length - 1 ] - _arr[0];
	},

	// variance of an array
	// flag indicates population vs sample
	variance : function( arr, flag ) {
		var mean = jStat.mean( arr ),
			stSum = 0,
			i = arr.length - 1;
		for ( ; i >= 0; i-- ) {
			stSum += Math.pow(( arr[i] - mean ), 2 );
		}
		return stSum / ( arr.length - ( flag ? 1 : 0 ));
	},

	// standard deviation of an array
	// flag indicates population vs sample
	stdev : function( arr, flag ) {
		return Math.sqrt( jStat.variance( arr, flag ));
	},

	// mean deviation (mean absolute deviation) of an array
	meandev : function( arr ) {
		var devSum = 0,
			mean = jStat.mean( arr ),
			i = arr.length - 1;
		for ( ; i >= 0; i-- ) {
			devSum += Math.abs( arr[i] - mean );
		}
		return devSum / arr.length;
	},


});



// exposing jStat
return jStat;

})( Math );
(function( jStat, Math ) {

// generate all distribution instance methods
(function( list ) {
	for ( var i = 0; i < list.length; i++ ) (function( func ) {
		// distribution instance method
		jStat[ func ] = function( a, b, c ) {
			if (!( this instanceof arguments.callee )) return new arguments.callee( a, b, c );
			this._a = a;
			this._b = b;
			this._c = c;
			return this;
		};
		// distribution method to be used on a jStat instance
		jStat.fn[ func ] = function( a, b, c ) {
			var newthis = jStat[ func ]( a, b, c );
			newthis.data = this;
			return newthis;
		};
		// sample instance method
		jStat[ func ].prototype.sample = function( arr ) {
			var a = this._a,
				b = this._b,
				c = this._c;
			if ( arr )
				return jStat.alter( arr, function() {
					return jStat[ func ].sample( a, b, c );
				});
			else
				return jStat[ func ].sample( a, b, c );
		};
		// generate the pdf, cdf and inv instance methods
		(function( vals ) {
			for ( var i = 0; i < vals.length; i++ ) (function( fnfunc ) {
				jStat[ func ].prototype[ fnfunc ] = function( x ) {
					var a = this._a,
						b = this._b,
						c = this._c;
					if ( isNaN( x )) {
						return jStat.fn.map.call( this.data, function( x ) {
							return jStat[ func ][ fnfunc ]( x, a, b, c );
						});
					}
					return jStat[ func ][ fnfunc ]( x, a, b, c );
				};
			})( vals[ i ]);
		})( 'pdf cdf inv'.split( ' ' ));
		// generate the mean, median, mode and variance instance methods
		(function( vals ) {
			for ( var i = 0; i < vals.length; i++ ) (function( fnfunc ) {
				jStat[ func ].prototype[ fnfunc ] = function() {
					return jStat[ func ][ fnfunc ]( this._a, this._b, this._c );
				};
			})( vals[ i ]);
		})( 'mean median mode variance'.split( ' ' ));
	})( list[ i ]);
})((
	'beta centralF cauchy chisquare exponential gamma invgamma kumaraswamy lognormal normal ' +
	'pareto studentt weibull uniform  binomial negbin hypgeom poisson triangular'
).split( ' ' ));


// extend beta function with static methods
jStat.extend( jStat.beta, {
	pdf : function( x, alpha, beta ) {
		return (x > 1 || x < 0) ? 0 : ( Math.pow( x, alpha - 1 ) * Math.pow( 1 - x, beta - 1 )) / jStat.betafn( alpha, beta );
	},

	cdf : function( x, alpha, beta ) {
		return (x > 1 || x < 0) ? (x > 1) * 1 : jStat.incompleteBeta( x, alpha, beta );
	},

	inv : function( x, alpha, beta ) {
		return jStat.incompleteBetaInv( x, alpha, beta );
	},

	mean : function( alpha, beta ) {
		return alpha / ( alpha + beta );
	},

	median : function( alpha, beta ) {
		// TODO: implement beta median
	},

	mode : function( alpha, beta ) {
		return ( alpha * beta ) / ( Math.pow( alpha + beta, 2 ) * ( alpha + beta + 1 ));
	},

	// return a random sample
	sample : function( alpha, beta ) {
		var u = jStat.randg( alpha );
		return u / ( u + jStat.randg( beta ));
	},

	variance : function( alpha, beta ) {
		return ( alpha * beta ) / ( Math.pow( alpha + beta, 2 ) * ( alpha + beta + 1 ));
	}
});





// extend gamma function with static methods
jStat.extend( jStat.gamma, {
	pdf : function( x, shape, scale ) {
		return Math.exp(( shape - 1 ) * Math.log( x ) - x / scale - jStat.gammaln( shape ) - shape * Math.log( scale ));
	},

	cdf : function( x, shape, scale ) {
		return jStat.gammap( shape, x / scale );
	},

	inv : function( p, shape, scale ) {
		return jStat.gammapinv( p, shape ) * scale;
	},

	mean : function( shape, scale ) {
		return shape * scale;
	},

	mode : function( shape, scale ) {
		if( shape > 1 ) return ( shape - 1 ) * scale;
		return undefined;
	},

	sample : function( shape, scale ) {
		return jStat.randg( shape ) * scale;
	},

	variance: function( shape, scale ) {
		return shape * scale * scale;
	}
});




// extend normal function with static methods
jStat.extend( jStat.normal, {
	pdf : function( x, mean, std ) {
		return Math.exp( -0.5 * Math.log( 2 * Math.PI ) - Math.log( std ) - Math.pow( x - mean, 2 ) / ( 2 * std * std ));
	},

	cdf : function( x, mean, std ) {
		return 0.5 * ( 1 + jStat.erf(( x - mean ) / Math.sqrt( 2 * std * std )));
	},

	inv : function( p, mean, std ) {
		return -1.41421356237309505 * std * jStat.erfcinv( 2 * p ) + mean;
	},

	mean : function( mean, std ) {
		return mean;
	},

	median : function( mean, std ) {
		return mean;
	},

	mode : function ( mean, std ) {
		return mean;
	},

	sample : function( mean, std ) {
		return jStat.randn() * std + mean;
	},

	variance : function( mean, std ) {
		return std * std;
	}
});



// extend uniform function with static methods
jStat.extend( jStat.binomial, {
	pdf : function( k, n, p ) {
		return ( p === 0 || p === 1 ) ?
			(( n * p ) === k ? 1 : 0 ) :
		jStat.combination( n, k ) * Math.pow( p, k ) * Math.pow( 1 - p, n - k );
	},

	cdf : function( x, n, p ) {
		var binomarr = [],
			k = 0;
		if ( x < 0 ) {
			return 0;
		}
		if ( x < n ) {
			for ( ; k <= x; k++ ) {
				binomarr[ k ] = jStat.binomial.pdf( k, n, p );
			}
			return jStat.sum( binomarr );
		}
		return 1;
	}
});



// extend negbin function with static methods
jStat.extend( jStat.negbin, {
	pdf : function( k, r, p ) {
		return k !== k | 0
			? false
		: k < 0
			? 0
		: jStat.combination( k + r - 1, k ) * Math.pow( 1 - p, r ) * Math.pow( p, k );
	},

	cdf : function( x, r, p ) {
		var sum = 0,
			k = 0;
		if ( x < 0 ) return 0;
		for ( ; k <= x; k++ ) {
			sum += jStat.negbin.pdf( k, r, p );
		}
		return sum;
	}
});




// extend uniform function with static methods
jStat.extend( jStat.poisson, {
	pdf : function( k, l ) {
		return Math.pow( l, k ) * Math.exp( -l ) / jStat.factorial( k );
	},

	cdf : function( x, l ) {
		var sumarr = [],
			k = 0;
		if ( x < 0 ) return 0;
		for ( ; k <= x; k++ ) {
			sumarr.push(jStat.poisson.pdf( k, l ));
		}
		return jStat.sum(sumarr);
	},

	mean : function( l ) {
		return l;
	},

	variance : function( l ) {
		return l;
	},

	sample : function( l ) {
		var p = 1, k = 0, L = Math.exp(-l);
		do {
			k++;
			p *= Math.random();
		} while (p > L);
		return k - 1;
	}
});


})( this.jStat, Math );
// Special functions //
(function( jStat, Math ) {


// extending static jStat methods
jStat.extend({

	// Log-gamma function
	gammaln : function( x ) {
		var j = 0,
			cof = [
				76.18009172947146, -86.50532032941677, 24.01409824083091,
				-1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5
			],
			ser = 1.000000000190015,
			xx, y, tmp;
		tmp = ( y = xx = x ) + 5.5;
		tmp -= ( xx + 0.5 ) * Math.log( tmp );
		for( ; j < 6; j++ ) ser += cof[j] / ++y;
		return Math.log( 2.5066282746310005 * ser / xx) - tmp;
	},

	// gamma of x
	gammafn : function( x ) {
		var p = [
				-1.716185138865495, 24.76565080557592, -379.80425647094563,
				629.3311553128184, 866.9662027904133, -31451.272968848367,
				-36144.413418691176, 66456.14382024054
			],
			q = [
				-30.8402300119739, 315.35062697960416, -1015.1563674902192,
				-3107.771671572311, 22538.118420980151, 4755.8462775278811,
				-134659.9598649693, -115132.2596755535
			],
			fact = false,
			n = 0,
			xden = 0,
			xnum = 0,
			y = x,
			i, z, yi, res, sum, ysq;
		if( y <= 0 ) {
			res = y % 1 + 3.6e-16;
			if ( res ) {
				fact = (!( y & 1 ) ? 1 : -1 ) * Math.PI / Math.sin( Math.PI * res );
				y = 1 - y;
			} else {
				return Infinity;
			}
		}
		yi = y;
		if ( y < 1 ) {
			z = y++;
		} else {
			z = ( y -= n = ( y | 0 ) - 1 ) - 1;
		}
		for ( i = 0; i < 8; ++i ) {
			xnum = ( xnum + p[i] ) * z;
			xden = xden * z + q[i];
		}
		res = xnum / xden + 1;
		if ( yi < y ) {
			res /= yi;
		} else if ( yi > y ) {
			for ( i = 0; i < n; ++i ) {
				res *= y;
				y++;
			}
		}
		if ( fact ) {
			res = fact / res;
		}
		return res;
	},

	// lower incomplete gamma function P(a,x)
	gammap : function( a, x ) {
		var aln = jStat.gammaln( a ),
			afn = jStat.gammafn( a ),
			ap = a,
			sum = 1 / a,
			del = sum,
			b = x + 1 - a,
			c = 1 / 1.0e-30,
			d = 1 / b,
			h = d,
			i = 1,
			afix = ( a >= 1 ) ? a : 1 / a,
			// calculate maximum number of itterations required for a
			ITMAX = -~( Math.log( afix ) * 8.5 + a * 0.4 + 17 ),
			an, endval;
		if ( x < 0 || a <= 0 ) {
			return NaN;
		} else if ( x < a + 1 ) {
			for ( ; i <= ITMAX; i++ ) {
				sum += del *= x / ++ap;
			}
			endval = sum * Math.exp( -x + a * Math.log( x ) - ( aln ));
		} else {
			for ( ; i <= ITMAX; i++ ) {
				an = -i * ( i - a );
				b += 2;
				d = an * d + b;
				c = b + an / c;
				d = 1 / d;
				h *= d * c;
			}
			endval = 1 - h * Math.exp( -x + a * Math.log( x ) - ( aln ));
		}
		return endval * afn / jStat.gammafn( a );
	},

	// natural log factorial of n
	factorialln : function( n ) {
		return n < 0 ? NaN : jStat.gammaln( n + 1 );
	},

	// factorial of n
	factorial : function( n ) {
		return n < 0 ? NaN : jStat.gammafn( n + 1 );
	},

	// combinations of n, m
	combination : function( n, m ) {
		// make sure n or m don't exceed the upper limit of usable values
		return ( n > 170 || m > 170 ) ?
			Math.exp( jStat.combinationln( n, m )) :
		( jStat.factorial( n ) / jStat.factorial( m )) / jStat.factorial( n - m );
	},

	combinationln : function( n, m ){
		return  jStat.factorialln( n ) - jStat.factorialln( m ) - jStat.factorialln( n - m );
	},

	// permutations of n, m
	permutation : function( n, m ) {
		return jStat.factorial( n ) / jStat.factorial( n - m );
	},

	// beta function
	betafn : function( x, y ) {
		// ensure arguments are positive
		if ( x <= 0 || y <= 0 ) return undefined;
		// make sure x + y doesn't exceed the upper limit of usable values
		return ( x + y > 170 ) ?
			Math.exp( jStat.betaln( x, y )) :
		jStat.gammafn( x ) * jStat.gammafn( y ) / jStat.gammafn( x + y );
	},

	// natural logarithm of beta function
	betaln : function( x, y ) {
		return jStat.gammaln( x ) + jStat.gammaln( y ) - jStat.gammaln( x + y );
	},

	// Returns the inverse incomplte gamma function
	gammapinv : function( p, a ) {
		var j = 0,
			a1 = a - 1,
			EPS = 1e-8,
			gln = jStat.gammaln( a ),
			x, err, t, u, pp, lna1, afac;
		if( p >= 1 ) return Math.max( 100, a + 100 * Math.sqrt( a ) );
		if( p <= 0 ) return 0;
		if( a > 1 ) {
			lna1 = Math.log( a1 );
			afac = Math.exp( a1 * ( lna1 - 1 ) - gln );
			pp = ( p < 0.5 ) ? p : 1 - p;
			t = Math.sqrt( -2 * Math.log( pp ));
			x = ( 2.30753 + t * 0.27061 ) / ( 1 + t * ( 0.99229 + t * 0.04481 )) - t;
			if( p < 0.5 ) x = -x;
			x = Math.max( 1e-3, a * Math.pow( 1 - 1 / ( 9 * a ) - x / ( 3 * Math.sqrt( a )), 3 ));
		} else {
			t = 1 - a * ( 0.253 + a * 0.12 );
			if( p < t ) x = Math.pow( p / t, 1 / a);
			else x = 1 - Math.log( 1 - ( p - t ) / ( 1 - t ));
		}
		for( ; j < 12; j++ ) {
			if( x <= 0 ) return 0;
			err = jStat.gammap( a, x ) - p;
			if( a > 1 ) t = afac * Math.exp( -( x - a1 ) + a1 * ( Math.log( x ) - lna1 ));
			else t = Math.exp( -x + a1 * Math.log( x ) - gln );
			u = err / t;
			x -= ( t = u / ( 1 - 0.5 * Math.min( 1, u * ( ( a - 1 ) / x - 1 ))));
			if( x <= 0 ) x = 0.5 * ( x + t );
			if( Math.abs( t ) < EPS * x ) break;
		}
		return x;
	},


});

// making use of static methods on the instance
(function( funcs ) {
	for ( var i = 0; i < funcs.length; i++ ) (function( passfunc ) {
		jStat.fn[ passfunc ] = function() {
			return jStat( jStat.map( this, function( value ) { return jStat[ passfunc ]( value ); }));
		};
	})( funcs[i] );
})( 'gammaln gammafn factorial factorialln'.split( ' ' ));

})( this.jStat, Math );

