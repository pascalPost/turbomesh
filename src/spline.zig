const std = @import("std");

pub fn FittingSpline(comptime dim: usize) type {
    return struct {
        const Self = @This();

        allocator: std.mem.Allocator,

        /// knots
        t: []f64,
        /// B-spline coefficients
        c: []f64,
        /// degree of spline
        k: u8,

        fn init(allocator: std.mem.Allocator, x: []const [dim]f64, k: u8) !Self {
            // set weights to one if no weights are passed in
            var w = try allocator.alloc(f64, x.len);
            @memset(w[0..], 1.0);
            defer allocator.free(w);

            if (dim >= 11) @compileError("fitting spline can only be up to 11 dimensions.");

            const iopt: c_int = 0; // new smoothing spline

            const ipar: c_int = 0; // arclength computed in function
            var u = try allocator.alloc(f64, x.len);
            @memset(u[0..], std.math.nan(f64));
            defer allocator.free(u);

            const idim: c_int = dim;
            const m: c_int = @intCast(x.len);
            const mx = m * idim;

            var ub: f64 = undefined;
            var ue: f64 = undefined;

            const k_int: c_int = k; // degree of spline
            const s: f64 = 0.0; // smoothing factor
            const nest: c_int = m + k_int + 1; // estimate number of knots
            var n: c_int = 0; // total number of knots

            var t = try allocator.alloc(f64, @intCast(nest)); // knots
            errdefer allocator.free(t);

            const nc = nest * idim;
            var c = try allocator.alloc(f64, @intCast(nc)); // coeffs
            errdefer allocator.free(c);

            var fp: f64 = undefined;

            const lwrk: c_int = m * (k_int + 1) + nest * (6 + idim + 3 * k_int);

            var wrk = try allocator.alloc(f64, @intCast(lwrk));
            defer allocator.free(wrk);

            // TODO perhaps we want to keep this buffer for later use (?)

            var iwrk = try allocator.alloc(c_int, @intCast(nest));
            defer allocator.free(iwrk);

            var ier: c_int = undefined;

            parcur(
                &iopt,
                &ipar,
                &idim,
                &m,
                u.ptr,
                &mx,
                &x[0],
                w.ptr,
                &ub,
                &ue,
                &k_int,
                &s,
                &nest,
                &n,
                t[0..].ptr,
                &nc,
                c[0..].ptr,
                &fp,
                wrk[0..].ptr,
                &lwrk,
                iwrk[0..].ptr,
                &ier,
            );

            if (ier > 0) {
                // TODO replace with logging
                std.debug.print("error in spline init (parcur); error code {}\n", .{ier});
                return error.SplineInit;
            }

            return Self{ .allocator = allocator, .t = t, .c = c, .k = k };
        }

        fn deinit(self: Self) void {
            self.allocator.free(self.t);
            self.allocator.free(self.c);
        }

        fn interpolate(self: *const Self, u: []const f64, values: [][dim]f64) !void {
            if (u.len != values.len) {
                std.debug.print("Mismatch of slice length ({} != {})", .{ u.len, values.len });
                return error.Mismatch;
            }

            const idim: c_int = dim;
            const n: c_int = @intCast(self.t.len);
            const nc: c_int = @intCast(self.c.len);
            const k_int: c_int = self.k;
            const m: c_int = @intCast(u.len);
            const mx: c_int = @as(c_int, @intCast(values.len)) * idim;
            var ier: c_int = undefined;
            curev(&idim, self.t.ptr, &n, self.c.ptr, &nc, &k_int, u.ptr, &m, &values[0], &mx, &ier);

            std.debug.assert(ier == 0);
        }
    };
}

// rename fortran underscored functions
const parcur = parcur_;
const curev = curev_;

/// given the ordered set of m points x(i) in the idim-dimensional space
/// and given also a corresponding set of strictly increasing values u(i)
/// and the set of positive numbers w(i),i=1,2,...,m, subroutine parcur
/// determines a smooth approximating spline curve s(u), i.e.
///     x1 = s1(u)
///     x2 = s2(u)       ub <= u <= ue
///     .........
///     xidim = sidim(u)
/// with sj(u),j=1,2,...,idim spline functions of degree k with common
/// knots t(j),j=1,2,...,n.
/// if ipar=1 the values ub,ue and u(i),i=1,2,...,m must be supplied by
/// the user. if ipar=0 these values are chosen automatically by parcur
/// as  v(1) = 0
///     v(i) = v(i-1) + dist(x(i),x(i-1)) ,i=2,3,...,m
///     u(i) = v(i)/v(m) ,i=1,2,...,m
///     ub = u(1) = 0, ue = u(m) = 1.
/// if iopt=-1 parcur calculates the weighted least-squares spline curve
/// according to a given set of knots.
/// if iopt>=0 the number of knots of the splines sj(u) and the position
/// t(j),j=1,2,...,n is chosen automatically by the routine. the smooth-
/// ness of s(u) is then achieved by minimalizing the discontinuity
/// jumps of the k-th derivative of s(u) at the knots t(j),j=k+2,k+3,...,
/// n-k-1. the amount of smoothness is determined by the condition that
/// f(p)=sum((w(i)*dist(x(i),s(u(i))))**2) be <= s, with s a given non-
/// negative constant, called the smoothing factor.
/// the fit s(u) is given in the b-spline representation and can be
/// evaluated by means of subroutine curev.
extern fn parcur_(
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
    u: [*]f64,
    // integer. on entry mx must specify the actual dimension of
    // the array x as declared in the calling (sub)program. mx must
    // not be too small (see x).
    mx: *const c_int,
    // real array of dimension at least idim*m.
    // before entry, x(idim*(i-1)+j) must contain the j-th coord-
    // inate of the i-th data point for i=1,2,...,m and j=1,2,...,
    // idim. unchanged on exit.
    x: [*]const f64,
    // real array of dimension at least (m). before entry, w(i)
    // must be set to the i-th value in the set of weights. the
    // w(i) must be strictly positive. unchanged on exit.
    // see also further comments.
    w: [*]const f64,
    // real values. on entry (in case ipar=1) ub and ue must
    // contain the lower and upper bound for the parameter u.
    // ub <=u(1), ue>= u(m). if ipar = 0 these values will
    // automatically be set to 0 and 1 by parcur.
    ub: *f64,
    ue: *f64,
    // integer. on entry k must specify the degree of the splines.
    // 1<=k<=5. it is recommended to use cubic splines (k=3).
    // the user is strongly dissuaded from choosing k even,together
    // with a small s-value. unchanged on exit.
    k: *const c_int,
    // real.on entry (in case iopt>=0) s must specify the smoothing
    // factor. s >=0. unchanged on exit.
    s: *const f64,
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
    n: *c_int,
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
    t: [*]f64,
    // integer. on entry nc must specify the actual dimension of
    // the array c as declared in the calling (sub)program. nc
    // must not be too small (see c). unchanged on exit.
    nc: *const c_int,
    // real array of dimension at least (nest*idim).
    // on succesful exit, this array will contain the coefficients
    // in the b-spline representation of the spline curve s(u),i.e.
    // the b-spline coefficients of the spline sj(u) will be given
    // in c(n*(j-1)+i),i=1,2,...,n-k-1 for j=1,2,...,idim.
    c: [*]f64,
    // real. unless ier = 10, fp contains the weighted sum of
    // squared residuals of the spline curve returned.
    fp: *f64,
    // real array of dimension at least m*(k+1)+nest*(6+idim+3*k).
    // used as working space. if the computation mode iopt=1 is
    // used, the values wrk(1),...,wrk(n) should be left unchanged
    // between subsequent calls.
    wrk: [*]f64,
    // integer. on entry,lwrk must specify the actual dimension of
    // the array wrk as declared in the calling (sub)program. lwrk
    // must not be too small (see wrk). unchanged on exit.
    lwrk: *const c_int,
    // integer array of dimension at least (nest).
    // used as working space. if the computation mode iopt=1 is
    // used,the values iwrk(1),...,iwrk(n) should be left unchanged
    // between subsequent calls.
    iwrk: [*]c_int,
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
    ier: *c_int,
) void;

/// curev evaluates in a number of points u(i),i=1,2,...,m
/// a spline curve s(u) of degree k and dimension idim, given in its
/// b-spline representation.
/// restrictions:
/// m >= 1
/// mx >= m*idim
/// t(k+1) <= u(i) <= u(i+1) <= t(n-k) , i=1,2,...,m-1.
extern fn curev_(
    /// integer, giving the dimension of the spline curve.
    idim: *const c_int,
    /// array,length n, which contains the position of the knots.
    t: [*]const f64,
    /// integer, giving the total number of knots of s(u).
    n: *const c_int,
    /// array,length nc, which contains the b-spline coefficients.
    c: [*]const f64,
    /// integer, giving the total number of coefficients of s(u).
    nc: *const c_int,
    /// integer, giving the degree of s(u).
    k: *const c_int,
    /// array,length m, which contains the points where s(u) must be evaluated.
    u: [*]const f64,
    /// integer, giving the number of points where s(u) must be evaluated.
    m: *const c_int,
    /// array,length mx,giving the value of s(u) at the different
    /// points. x(idim*(i-1)+j) will contain the j-th coordinate
    /// of the i-th point on the curve.
    x: [*]f64,
    /// integer, giving the dimension of the array x. mx >= m*idim
    mx: *const c_int,
    /// error flag
    /// ier = 0 : normal return
    /// ier =10 : invalid input data (see restrictions)
    ier: *c_int,
) void;

test "spline" {
    const allocator = std.testing.allocator;
    const dim = 2;
    const order = 3;

    const x = [_][dim]f64{ .{ 0.0, 0.0 }, .{ 0.5, 0.5 }, .{ 1.0, 1.0 }, .{ 2.0, 2.0 }, .{ 3.0, 3.0 }, .{ 4.0, 4.0 } };

    const spline = try FittingSpline(dim).init(allocator, x[0..], order);
    defer spline.deinit();

    {
        const u = [_]f64{ 0.0, 0.125, 0.25, 0.5, 0.75, 1.0 };
        var values: [u.len][dim]f64 = undefined;
        try spline.interpolate(&u, &values);
        for (values, x) |v_i, x_i| {
            for (v_i, x_i) |v_ij, x_ij| try std.testing.expectApproxEqAbs(x_ij, v_ij, 1e-15);
        }
    }
}
