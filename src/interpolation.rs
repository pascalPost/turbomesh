// Copyright (c) 2022 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

use dierckx_sys::curev_;
use libc::{c_double, c_int};

extern "C" {
    // given the ordered set of m points x(i) in the idim-dimensional space
    // and given also a corresponding set of strictly increasing values u(i)
    // and the set of positive numbers w(i),i=1,2,...,m, subroutine parcur
    // determines a smooth approximating spline curve s(u), i.e.
    //     x1 = s1(u)
    //     x2 = s2(u)       ub <= u <= ue
    //     .........
    //     xidim = sidim(u)
    // with sj(u),j=1,2,...,idim spline functions of degree k with common
    // knots t(j),j=1,2,...,n.
    // if ipar=1 the values ub,ue and u(i),i=1,2,...,m must be supplied by
    // the user. if ipar=0 these values are chosen automatically by parcur
    // as  v(1) = 0
    //     v(i) = v(i-1) + dist(x(i),x(i-1)) ,i=2,3,...,m
    //     u(i) = v(i)/v(m) ,i=1,2,...,m
    //     ub = u(1) = 0, ue = u(m) = 1.
    // if iopt=-1 parcur calculates the weighted least-squares spline curve
    // according to a given set of knots.
    // if iopt>=0 the number of knots of the splines sj(u) and the position
    // t(j),j=1,2,...,n is chosen automatically by the routine. the smooth-
    // ness of s(u) is then achieved by minimalizing the discontinuity
    // jumps of the k-th derivative of s(u) at the knots t(j),j=k+2,k+3,...,
    // n-k-1. the amount of smoothness is determined by the condition that
    // f(p)=sum((w(i)*dist(x(i),s(u(i))))**2) be <= s, with s a given non-
    // negative constant, called the smoothing factor.
    // the fit s(u) is given in the b-spline representation and can be
    // evaluated by means of subroutine curev.
    pub fn parcur_(
        // integer flag. on entry iopt must specify whether a weighted
        // least-squares spline curve (iopt=-1) or a smoothing spline
        // curve (iopt=0 or 1) must be determined.if iopt=0 the routine
        // will start with an initial set of knots t(i)=ub,t(i+k+1)=ue,
        // i=1,2,...,k+1. if iopt=1 the routine will continue with the
        // knots found at the last call of the routine.
        // attention: a call with iopt=1 must always be immediately
        // preceded by another call with iopt=1 or iopt=0. unchanged on exit.
        iopt: *const c_int,
        // integer flag. on entry ipar must specify whether (ipar=1)
        // the user will supply the parameter values u(i), ub and ue
        // or whether (ipar=0) these values are to be calculated by
        // parcur. unchanged on exit.
        ipar: *const c_int,
        // integer. on entry idim must specify the dimension of the
        // curve. 0 < idim < 11. unchanged on exit.
        idim: *const c_int,
        // integer. on entry m must specify the number of data points.
        // m > k. unchanged on exit.
        m: *const c_int,
        // real array of dimension at least (m). in case ipar=1, before
        // entry, u(i) must be set to the i-th value of the parameter
        // variable u for i=1,2,...,m. These values must then be
        // supplied in strictly ascending order and will be unchanged
        // on exit. in case ipar=0, on exit, array u will contain the
        // values u(i) as determined by parcur.
        u: *mut c_double,
        // integer. on entry mx must specify the actual dimension of
        // the array x as declared in the calling (sub)program. mx must
        // not be too small (see x).
        mx: *const c_int,
        // real array of dimension at least idim*m.
        // before entry, x(idim*(i-1)+j) must contain the j-th coord-
        // inate of the i-th data point for i=1,2,...,m and j=1,2,...,
        // idim. unchanged on exit.
        x: *const c_double,
        // real array of dimension at least (m). before entry, w(i)
        // must be set to the i-th value in the set of weights. the
        // w(i) must be strictly positive. unchanged on exit.
        // see also further comments.
        w: *const c_double,
        // real values. on entry (in case ipar=1) ub and ue must
        // contain the lower and upper bound for the parameter u.
        // ub <=u(1), ue>= u(m). if ipar = 0 these values will
        // automatically be set to 0 and 1 by parcur.
        ub: *mut c_double,
        ue: *mut c_double,
        // integer. on entry k must specify the degree of the splines.
        // 1<=k<=5. it is recommended to use cubic splines (k=3).
        // the user is strongly dissuaded from choosing k even,together
        // with a small s-value. unchanged on exit.
        k: *const c_int,
        // real.on entry (in case iopt>=0) s must specify the smoothing
        // factor. s >=0. unchanged on exit.
        s: *const c_double,
        // integer. on entry nest must contain an over-estimate of the
        // total number of knots of the splines returned, to indicate
        // the storage space available to the routine. nest >=2*k+2.
        // in most practical situation nest=m/2 will be sufficient.
        // always large enough is nest=m+k+1, the number of knots
        // needed for interpolation (s=0). unchanged on exit.
        nest: *const c_int,
        // integer.
        // unless ier = 10 (in case iopt >=0), n will contain the
        // total number of knots of the smoothing spline curve returned
        // if the computation mode iopt=1 is used this value of n
        // should be left unchanged between subsequent calls.
        // in case iopt=-1, the value of n must be specified on entry.
        n: *mut c_int,
        // real array of dimension at least (nest).
        // on succesful exit, this array will contain the knots of the
        // spline curve, i.e. the position of the interior knots t(k+2),
        // t(k+3),..,t(n-k-1) as well as the position of the additional
        // t(1)=t(2)=...=t(k+1)=ub and t(n-k)=...=t(n)=ue needed for
        // the b-spline representation.
        // if the computation mode iopt=1 is used, the values of t(1),
        // t(2),...,t(n) should be left unchanged between subsequent
        // calls. if the computation mode iopt=-1 is used, the values
        // t(k+2),...,t(n-k-1) must be supplied by the user, before
        // entry. see also the restrictions (ier=10).
        t: *mut c_double,
        // integer. on entry nc must specify the actual dimension of
        // the array c as declared in the calling (sub)program. nc
        // must not be too small (see c). unchanged on exit.
        nc: *const c_int,
        // real array of dimension at least (nest*idim).
        // on succesful exit, this array will contain the coefficients
        // in the b-spline representation of the spline curve s(u),i.e.
        // the b-spline coefficients of the spline sj(u) will be given
        // in c(n*(j-1)+i),i=1,2,...,n-k-1 for j=1,2,...,idim.
        c: *mut c_double,
        // real. unless ier = 10, fp contains the weighted sum of
        // squared residuals of the spline curve returned.
        fp: *mut c_double,
        // real array of dimension at least m*(k+1)+nest*(6+idim+3*k).
        // used as working space. if the computation mode iopt=1 is
        // used, the values wrk(1),...,wrk(n) should be left unchanged
        // between subsequent calls.
        wrk: *mut c_double,
        // integer. on entry,lwrk must specify the actual dimension of
        // the array wrk as declared in the calling (sub)program. lwrk
        // must not be too small (see wrk). unchanged on exit.
        lwrk: *const c_int,
        // integer array of dimension at least (nest).
        // used as working space. if the computation mode iopt=1 is
        // used,the values iwrk(1),...,iwrk(n) should be left unchanged
        // between subsequent calls.
        iwrk: *mut c_int,
        // integer. unless the routine detects an error, ier contains a
        // non-positive value on exit, i.e.
        //  ier=0  : normal return. the curve returned has a residual sum of
        //           squares fp such that abs(fp-s)/s <= tol with tol a relat-
        //           ive tolerance set to 0.001 by the program.
        //  ier=-1 : normal return. the curve returned is an interpolating
        //           spline curve (fp=0).
        //  ier=-2 : normal return. the curve returned is the weighted least-
        //           squares polynomial curve of degree k.in this extreme case
        //           fp gives the upper bound fp0 for the smoothing factor s.
        //  ier=1  : error. the required storage space exceeds the available
        //           storage space, as specified by the parameter nest.
        //           probably causes : nest too small. if nest is already
        //           large (say nest > m/2), it may also indicate that s is
        //           too small
        //           the approximation returned is the least-squares spline
        //           curve according to the knots t(1),t(2),...,t(n). (n=nest)
        //           the parameter fp gives the corresponding weighted sum of
        //           squared residuals (fp>s).
        //  ier=2  : error. a theoretically impossible result was found during
        //           the iteration proces for finding a smoothing spline curve
        //           with fp = s. probably causes : s too small.
        //           there is an approximation returned but the corresponding
        //           weighted sum of squared residuals does not satisfy the
        //           condition abs(fp-s)/s < tol.
        //  ier=3  : error. the maximal number of iterations maxit (set to 20
        //           by the program) allowed for finding a smoothing curve
        //           with fp=s has been reached. probably causes : s too small
        //           there is an approximation returned but the corresponding
        //           weighted sum of squared residuals does not satisfy the
        //           condition abs(fp-s)/s < tol.
        //  ier=10 : error. on entry, the input data are controlled on validity
        //           the following restrictions must be satisfied.
        //           -1<=iopt<=1, 1<=k<=5, m>k, nest>2*k+2, w(i)>0,i=1,2,...,m
        //           0<=ipar<=1, 0<idim<=10, lwrk>=(k+1)*m+nest*(6+idim+3*k),
        //           nc>=nest*idim
        //           if ipar=0: sum j=1,idim (x(idim*i+j)-x(idim*(i-1)+j))**2>0
        //                      i=1,2,...,m-1.
        //           if ipar=1: ub<=u(1)<u(2)<...<u(m)<=ue
        //           if iopt=-1: 2*k+2<=n<=min(nest,m+k+1)
        //                       ub<t(k+2)<t(k+3)<...<t(n-k-1)<ue
        //                          (ub=0 and ue=1 in case ipar=0)
        //                     the schoenberg-whitney conditions, i.e. there
        //                     must be a subset of data points uu(j) such that
        //                       t(j) < uu(j) < t(j+k+1), j=1,2,...,n-k-1
        //           if iopt>=0: s>=0
        //                       if s=0 : nest >= m+k+1
        //           if one of these conditions is found to be violated,control
        //           is immediately repassed to the calling program. in that
        //           case there is no approximation returned.
        ier: *mut c_int,
    );
}

pub struct FittingSpline<const DIM: usize> {
    /// knots
    t: Box<[f64]>,
    /// B-spline coefficients
    c: Box<[f64]>,
    /// degree of spline
    k: u8,
}

impl<const DIM: usize> FittingSpline<DIM> {
    pub fn new(
        x: &[[f64; DIM]],
        w_opt: Option<&[f64]>,
        u_opt: Option<(&[f64], f64, f64)>,
        k: u8,
    ) -> Result<Self, i32> {
        let res = Self::prepare(x, w_opt, u_opt, k);

        match res {
            Ok((data, _, _, _)) => Ok(data),
            Err(e) => Err(e),
        }
    }

    fn prepare(
        x: &[[f64; DIM]],
        w_opt: Option<&[f64]>,
        u_opt: Option<(&[f64], f64, f64)>,
        k: u8,
    ) -> Result<(Self, Option<Box<[f64]>>, f64, i32), i32> {
        let n = x.len();

        // set weights to one if no weights are passed in
        let w_mem: Option<Vec<f64>>;
        let w = if let Some(w) = w_opt {
            w
        } else {
            w_mem = Some(vec![1.0; n]);
            &w_mem.as_ref().unwrap()[..]
        };

        let mut u_mem: Option<Vec<f64>> = None;
        let (ipar, u, mut ub, mut ue) = if let Some((u, ub, ue)) = u_opt {
            // ipar=1 : user suplied u, ub, ue
            assert!(u.len() == n);
            assert!(ub <= u[0]);
            assert!(ue >= u[n - 1]);
            let ipar = 1;
            (ipar, u, ub, ue)
        } else {
            // ipar=0 : u computed in the function
            u_mem = Some(vec![f64::NAN; n]);
            let ipar = 0;
            (ipar, &u_mem.as_ref().unwrap()[..], 0.0, 1.0)
        };

        unsafe {
            let iopt: c_int = 0;
            let idim = DIM as c_int;

            assert!(idim > 0 && idim < 11, "0 < idim < 11 must hold");

            let m: c_int = n as c_int; // number of data points

            let mx: c_int = m * idim;

            let s: c_double = 0.0;
            let nest: c_int = m + k as c_int + 1;
            let mut n: c_int = 0;

            let mut t = vec![0.0; nest as usize];

            let nc: c_int = nest * idim;
            let mut c = vec![0.0; nc as usize];

            let mut fp: c_double = 0.0;

            let lwrk: c_int = m * (k as c_int + 1) + nest * (6 + idim + 3 * k as c_int);
            let mut wrk = vec![0.0; lwrk as usize];
            let mut iwrk = vec![0; nest as usize];

            let mut ier: c_int = 0;

            parcur_(
                &iopt,
                &ipar,
                &idim,
                &m,
                u.as_ptr() as *mut c_double,
                &mx,
                x.as_ptr() as *const c_double,
                w.as_ptr(),
                &mut ub,
                &mut ue,
                &(k as c_int),
                &s,
                &nest,
                &mut n,
                t.as_mut_ptr(),
                &nc,
                c.as_mut_ptr(),
                &mut fp,
                wrk.as_mut_ptr(),
                &lwrk,
                iwrk.as_mut_ptr(),
                &mut ier,
            );

            // assert!(ier <= 0, "no error, but ier == {} != 0", ier);

            // assert!(ub == 0.0);
            // assert!(ue == 1.0);

            if ier <= 0 {
                Ok((
                    Self {
                        t: t.into_boxed_slice(),
                        c: c.into_boxed_slice(),
                        k,
                    },
                    u_mem.map(|v| v.into_boxed_slice()),
                    fp,
                    ier,
                ))
            } else {
                Err(ier)
            }
        }
    }

    pub fn interpolate(&self, u: &[f64]) -> Box<[[f64; DIM]]> {
        let mut res = vec![[f64::NAN; DIM]; u.len()];

        unsafe {
            let idim = DIM as c_int;
            let mut ier: c_int = 0;
            curev_(
                &idim,
                self.t.as_ptr(),
                &(self.t.len() as c_int),
                self.c.as_ptr(),
                &(self.c.len() as c_int),
                &(self.k as c_int),
                u.as_ptr() as *const c_double,
                &(u.len() as c_int),
                res.as_mut_ptr() as *mut c_double,
                &((res.len() * DIM as usize) as c_int),
                &mut ier,
            );

            assert!(ier == 0);
        }

        res.into_boxed_slice()
    }
}
