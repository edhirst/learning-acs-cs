"""F^-1 on uniformly distributed points of the round S^4 x T^2: residual |F(F^-1(q)) - q|, the initial-guess quality
(the explicit chain is exact outside the linear-blend shells), and the transported J evaluated there.
Run from the repo root: python3 extra/s4t2_explicit/check_inverse.py [n]"""
from __future__ import annotations
import math, os, sys, time
import torch
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
torch.set_default_dtype(torch.float64)
from extra.s4t2_explicit import tinverse as TI, nijenhuis_round as NR


def check_transport_base(algebraic_radial=False,algebraic_profile=False,finite_height=False):
    """Independent unfitted base: certified shell roots, axes and round targets."""
    from extra.s4t2_explicit import tcore as T, tlift as L
    from torch.func import jvp
    gen = torch.Generator().manual_seed(623)
    targets = []
    sphere = torch.randn(2000,5,generator=gen)
    sphere /= sphere.norm(dim=-1,keepdim=True)
    z1,z2 = torch.complex(sphere[:,0],sphere[:,1]),torch.complex(sphere[:,2],sphere[:,3])
    hh = 2*z1*z2.conj()
    targets.append((torch.stack([hh.real,hh.imag,T.abs2(z1)-T.abs2(z2)],-1),sphere[:,4]))
    worst_shell = 0.
    C0,Ci = T.height_cusps(finite_height)
    for C,infinity in ((C0,False),(Ci,True)):
        # The forward base's second derivatives at the node must agree with
        # its quadratic Hopf model, not differentiate a norm at zero.
        origin = torch.zeros(1,4)
        basis = torch.eye(4)
        def at_node(z):
            vv,hh = L.transport_base_chart(C,torch.complex(z[:,0],z[:,1]),torch.complex(z[:,2],z[:,3]),algebraic_radial=algebraic_radial,algebraic_profile=algebraic_profile)
            return torch.cat([vv,hh[:,None]],-1)
        def model(z):
            scale = L.SINF if infinity else L.S0
            xx,xi,yy,yi = z.unbind(-1)
            pp = torch.stack([2*math.sqrt(C.a0)*(xx*yy-xi*yi),
                              -2*math.sqrt(C.a0)*(xx*yi+xi*yy),
                              C.a0*(xx*xx+xi*xi)-yy*yy-yi*yi],-1)
            vv = pp@(L.MINF if infinity else L.M0).T
            rad = scale*(C.a0*(xx*xx+xi*xi)+yy*yy+yi*yi)
            if algebraic_radial:
                vv = 2*vv/(1+rad*rad)[:,None]
                hh = (1 if infinity else -1)*(1-rad)/torch.sqrt(1+rad*rad)
            else:
                hh = (1 if infinity else -1)*torch.sqrt(1-rad)
            return torch.cat([vv,hh[:,None]],-1)
        for i in range(4):
            for j in range(4):
                second = lambda f: jvp(lambda z:jvp(f,(z,),(basis[i:i+1],))[1],(origin,),(basis[j:j+1],))[1]
                actual,expected = second(at_node),second(model)
                assert torch.isfinite(actual).all() and torch.allclose(actual,expected,atol=1e-13,rtol=1e-13)
        pp = torch.randn(1024,3,generator=gen)
        pp /= pp.norm(dim=-1,keepdim=True)
        rr = torch.linspace(0,L.EPS2,1024)
        pp *= rr[:,None]
        vv = L.transport_base_local(pp,infinity,algebraic_profile=algebraic_profile)
        inv,inside = TI.transport_shell_inverse(vv,infinity,algebraic_profile=algebraic_profile)
        err = (inv-pp).abs().max().item()
        assert err < 1e-12 and inside[:-1].all()
        worst_shell = max(worst_shell,err)
        nv = vv.norm(dim=-1)
        if algebraic_radial:
            targets.append((2*vv/(1+nv*nv)[:,None],(1 if infinity else -1)*(1-nv)/torch.sqrt(1+nv*nv)))
        else:
            targets.append((vv,(1 if infinity else -1)*torch.sqrt(1-nv)))
        # Central axes, including the height half-cell and the node-blend endpoints.
        radii = torch.cat([torch.tensor([0.,.004,.008,.02,.05,1.]),torch.linspace(.001,.999,122)])
        xx = torch.cat([torch.sqrt(radii/C.a0),torch.zeros_like(radii)]).to(torch.complex128)
        yy = torch.cat([torch.zeros_like(radii),torch.sqrt(radii)]).to(torch.complex128)
        targets.append(L.transport_base_chart(C,xx,yy,algebraic_radial=algebraic_radial,algebraic_profile=algebraic_profile))
    v,X = (torch.cat([pair[k] for pair in targets]) for k in range(2))
    result = TI.transport_base_inverse(v,X,algebraic_radial=algebraic_radial,algebraic_profile=algebraic_profile,finite_height=finite_height)
    assert torch.equal(torch.cat([idx for idx,_ in result.values()]).sort().values,torch.arange(len(X)))
    error = 0.
    for key,(idx,data) in result.items():
        if key == 'mid':
            t = torch.complex(data[:,0],data[:,1])
            vv,xx = L.transport_base_middle(t,data[:,2],data[:,3],torch.log(t),algebraic_radial=algebraic_radial,algebraic_profile=algebraic_profile,finite_height=finite_height)
        else:
            vv,xx = L.transport_base_chart(C0 if key == 'c0' else Ci,
                    torch.complex(data[:,0],data[:,1]),torch.complex(data[:,2],data[:,3]),algebraic_radial=algebraic_radial,algebraic_profile=algebraic_profile)
        assert torch.isfinite(data).all() and torch.isfinite(vv).all() and torch.isfinite(xx).all()
        error = max(error,(vv-v[idx]).abs().max().item(),(xx-X[idx]).abs().max().item())
    assert error < 1e-10
    print(f'Transport base (algebraic radial={algebraic_radial}, profile={algebraic_profile}, finite height={finite_height}): {len(X)} random/shell/pole/axis targets; shell error {worst_shell:.1e}, full error {error:.1e}: PASS')
    print('Transport base: all 16 second-derivative directions at each exact pole: PASS')


def check_algebraic_axis_jets():
    """All derivatives through order three at tc=0 vs independent Taylor polynomials."""
    from extra.s4t2_explicit import tlift as L
    from torch.func import jvp
    import itertools
    A,B = (4*math.pi)**2,math.exp(4)
    origin = torch.zeros(1,2)
    directions = torch.eye(2)
    for infinity in (False,True):
        phase = torch.exp(torch.tensor(2j*math.pi*((-.17+.5) if infinity else .17)))
        def actual(z):
            u = L.algebraic_hopf_cusp(torch.complex(z[:,0],z[:,1]),torch.full((len(z),),.17),infinity)
            return torch.view_as_real(u).flatten(1)
        def taylor(z):
            t = torch.complex(z[:,0],z[:,1])
            r2 = (z*z).sum(-1)
            if infinity:
                u1 = phase*(1-A/B*r2/2)
                u2 = phase*math.sqrt(A/B)*t*(1+(1/A-2*A/B)*r2/2)
            else:
                u1 = phase*math.sqrt(A)*t*(1+(B/A-2*A)*r2/2)
                u2 = phase*(1-A*r2/2)
            return torch.view_as_real(torch.stack([u1,u2],-1)).flatten(1)
        worst = 0.
        for degree in range(1,4):
            for inds in itertools.product(range(2),repeat=degree):
                def derivative(fun):
                    for idx in inds:
                        previous = fun
                        fun = lambda z,previous=previous,idx=idx:jvp(previous,(z,),(directions[idx:idx+1],))[1]
                    return fun(origin)
                err = (derivative(actual)-derivative(taylor)).abs().max().item()
                worst = max(worst,err)
                assert err < 1e-10
        print(f'Algebraic profile cusp {"inf" if infinity else "0"}: 14 central-axis derivative checks through order 3; error {worst:.1e}: PASS')

def check_elementary_closed_blocks():
    """Independent numerical checks of programme 11; NOT the huge-count map.

    Run alone with --closed. Full-count nonvanishing is the written certificate;
    the practical counts below check formulas and local overlaps only.
    """
    import mpmath as mp
    import numpy as np
    from extra.s4t2_explicit import cusp as C, tcore as T
    with mp.workdps(55):
        worst = mp.mpf(0)
        count = 0
        phi = (mp.mpf('.17'),mp.mpf('-.23'))
        for inf in (False,True):
            nu = mp.exp(-2) if inf else mp.mpf(1)
            for r in map(mp.mpf,('.001','.1','.299','.3','.32','.4','.49')):
                for a in map(mp.mpf,('-.499','-.2','0','.2','.499')):
                    tc = r*mp.exp(mp.j*mp.mpf('.7'))
                    rho = C.elementary_cusp_logradius(a,tc,nu,dps=55)
                    y = mp.exp(rho)
                    x = tc/y
                    t = 1/tc if inf else tc
                    lt = mp.log(t)
                    ss = lt/(2*mp.j*mp.pi)
                    zm,wm = (-mp.j*rho,0) if inf else (-mp.j*rho+mp.pi*ss*(ss-1),2*t*ss)
                    lattice = C.scalar_lattice_coords(zm,wm,t,lt)
                    data = (t,lattice[0],lattice[3],lt,a)
                    mid,norm = C.elliptic_bundle_embedding('mid',data,dps=55,with_norm=True,elementary=True)
                    # Drop the evaluator's guard digits before comparing an
                    # exactly attained endpoint (beta=1/2, norm^2=1/2).
                    norm = +norm
                    assert mp.mpf('.5')<=norm<=2
                    if r<=mp.mpf('.3'):
                        chart = C.elliptic_bundle_embedding('ci' if inf else 'c0',(x,y,a),dps=55,elementary=True)
                        worst = max(worst,mp.norm(mp.matrix(mid)-mp.matrix(chart)))
                    # Branch change and circle action, using an independent point.
                    shifted = (t,lattice[0]-lattice[3],lattice[3],lt+2*mp.pi*mp.j,a)
                    other = C.elliptic_bundle_embedding('mid',shifted,dps=55,elementary=True)
                    worst = max(worst,mp.norm(mp.matrix(mid)-mp.matrix(other)))
                    phase = mp.exp(mp.j*mp.mpf('.37'))
                    moved = (t,lattice[0]+mp.mpf('.37')/(2*mp.pi),lattice[3],lt,a)
                    other = C.elliptic_bundle_embedding('mid',moved,dps=55,elementary=True)
                    worst = max(worst,mp.norm(mp.matrix(other)-phase*mp.matrix(mid)))
                    deck = (t,lattice[0],lattice[3]+1,lt,a+(-1 if inf else 1))
                    other = C.elliptic_bundle_embedding('mid',deck,dps=55,elementary=True)
                    worst = max(worst,mp.norm(mp.matrix(other)-mp.matrix(mid)))
                    # Restore T2, then undo it by independent torch lattice/slice
                    # formulas and the previously tested finite-height class.
                    kind,h = C.elementary_restore_torus('mid',data,phi)
                    zt,wt,tt = [torch.tensor([complex(v)]) for v in h]
                    ltt = torch.tensor([complex(lt)])
                    lx = T.lattice_coords(zt,wt,tt,ltt)[0]
                    cc = T.FINITE_CUSPINF if inf else T.FINITE_CUSP0
                    convert = T.to_cuspinf if inf else T.to_cusp0
                    zz,ss_t = convert(tt,lx[0:1],lx[3:4],ltt)
                    HH = cc.mu_torus(zz,ss_t)[0]
                    ct = torch.stack(cc.c_tilde(ss_t),-1)[0]
                    if inf:
                        correction = (HH+lx[3])*ct
                    else:
                        xi23 = torch.stack(T.xi23(tt),-1)[0]
                        xi4 = T.xi4(tt)[0]
                        correction = xi23+(HH-lx[3]+xi4)*ct
                    beta = 1-T.bump_step((ss_t.abs()-.3)/.2)[0]
                    recovered = lx[1:3]-beta*correction
                    assert (recovered-torch.tensor(list(map(float,phi)))).abs().max()<2e-10
                    count += 1
        assert worst<mp.mpf('1e-40')
        print(f'Closed blocks: {count} cusp/collar embedding, branch, circle and independent torus checks; error {float(worst):.1e}: PASS',flush=True)
        # The terminal-cut moment has an elementary inverse on every height
        # cell, including the collar and the two flat transition seams.
        height_error = mp.mpf(0)
        height_cases = 0
        for inf in (False,True):
            nu = mp.exp(-2) if inf else mp.mpf(1)
            for radius in map(mp.mpf,('1e-12','.3','.31','.4','.49','.5')):
                tc = radius if radius==mp.mpf('.5') else radius*mp.exp(mp.j*mp.mpf('.7'))
                for a in map(mp.mpf,('-.5','-.4','-.2','-.1','0','.1','.2','.4','.5','1.2')):
                    rho = C.elementary_cusp_logradius(a,tc,nu,dps=55)
                    recovered = C.elementary_cusp_height(rho,tc,nu,dps=55)
                    height_error = max(height_error,abs(recovered-a))
                    assert abs(recovered-a)<mp.mpf('1e-40')
                    height_cases += 1
                    if radius<mp.mpf('.5'):
                        y = mp.exp(rho)
                        chart_height = C.elementary_cusp_height_chart(tc/y,y,nu,dps=55)
                        height_error = max(height_error,abs(chart_height-a))
                        assert abs(chart_height-a)<mp.mpf('1e-40')
                        height_cases += 1
            P = mp.log(nu*mp.mpf('.4'))-2*mp.mpf('.4')**2/(1-mp.mpf('.4')**2)
            q = mp.exp(2*P)
            lam = C.scalar_flat_step(mp.mpf('.5'))
            for seam in map(mp.mpf,('.1','.2','.4')):
                for gap in map(mp.mpf,('1e-6','1e-10','1e-40')):
                    for sign in (-1,1):
                        a = seam+sign*gap
                        value = C.inverse_height_moment(a,q,lam)
                        recovered = C.elementary_forward_moment(value,q,lam)
                        height_error = max(height_error,abs(recovered-a))
                        assert abs(recovered-a)<mp.mpf('1e-40')
                        height_cases += 1
            for a in map(mp.mpf,('-.55','-.5','-.49','-.2','.2','.49','.5','.55')):
                x,y = C.elementary_cusp_reference(0,a,nu)
                recovered = C.elementary_cusp_height_chart(x,y,nu,dps=55)
                height_error = max(height_error,abs(recovered-a))
                assert abs(recovered-a)<mp.mpf('1e-40')
                height_cases += 1
            assert C.elementary_cusp_height_chart(0,0,nu,dps=55)==0
            height_cases += 1
            assert C.elementary_cusp_height_chart(mp.mpf('.5'),mp.mpf(1),nu,dps=55)==-mp.mpf('.5')
            height_cases += 1
        print(f'Elementary two-way height: {height_cases} cell, seam, axis and node '
              f'compositions; error {float(height_error):.1e}: PASS',flush=True)
        # Reverse the complete closed map, including tau_f, both pole caps,
        # active cusp collars, the middle phase, and source chart changes.
        forward_error = mp.mpf(0)
        forward_cases = 0
        phi = (mp.mpf('.17'),mp.mpf('.29'))
        xi = (mp.sqrt(mp.mpf('.7')),
              mp.sqrt(mp.mpf('.3'))*mp.exp(mp.j*mp.mpf('.7')))
        def compare_forward(kind,data,Z,X,angles,steps):
            nonlocal forward_error,forward_cases
            recovered,rx,ra = C.elementary_global_forward(
                kind,data,projection_steps=steps,phase_steps=32,dps=55)
            error = max(abs(rx-X),*(abs(recovered[j]-Z[j]) for j in range(2)),
                        *(abs(ra[j]-angles[j]) for j in range(2)))
            forward_error = max(forward_error,error)
            forward_cases += 1
            assert error<mp.mpf('1e-40')
        for X in map(mp.mpf,('-1','-.8','-.15','0','.15','.8','1')):
            Z = tuple(mp.sqrt(1-X*X)*v for v in xi)
            kind,data = C.elementary_global_coordinates(
                Z,X,phi,projection_steps=64,phase_steps=32,dps=55)
            compare_forward(kind,data,Z,X,phi,64)
        for sign in (-1,1):
            for gap,steps in ((mp.mpf('4e-6'),64),(mp.mpf('6e-6'),512)):
                X = sign*(1-gap)
                Z = tuple(mp.sqrt(1-X*X)*v for v in xi)
                kind,data = C.elementary_global_coordinates(
                    Z,X,phi,projection_steps=steps,phase_steps=32,dps=55)
                compare_forward(kind,data,Z,X,phi,steps)
                if gap==mp.mpf('4e-6'):
                    # The same cap point also has a middle holomorphic chart.
                    x,y,w = data
                    tc = x*y
                    if kind=='c0':
                        t = tc
                        s = mp.log(t)/(2*mp.j*mp.pi)
                        zm = -mp.j*mp.log(y)+mp.pi*s*(s-1)
                        wm = w+2*t*s
                    else:
                        assert kind=='ci'
                        t = 1/tc
                        wm = w*t*t
                        zm = -mp.j*mp.log(y)+mp.j*wm/t
                    compare_forward('mid',(zm,wm,t),Z,X,phi,steps)
        # A cusp W point has deck representatives with H_f shifted by n.
        for inf in (False,True):
            nu = mp.exp(-2) if inf else mp.mpf(1)
            X = mp.mpf('.8') if inf else mp.mpf('-.8')
            Z = tuple(mp.sqrt(1-X*X)*v for v in xi)
            kind,reference = C.elementary_global_slice(
                Z,X,projection_steps=64,phase_steps=32,dps=55)
            xs,ys = reference[:2]
            t = xs*ys
            h = nu*mp.exp(-2*abs(t)**2/(1-abs(t)**2))*t
            e2,e3,e4,a2,a3 = C._elementary_cusp_basis(t,inf)
            c = C._elementary_real_coeffs(e4,e2,e3)
            H0 = C.scalar_finite_height(xs,ys,nu)
            for deck in (-1,1):
                xn,yn = xs*h**(-deck),ys*h**deck
                coeff = tuple(phi[j]+(H0+deck)*c[j] for j in range(2))
                A = coeff[0]*a2+coeff[1]*a3
                source = (xn*mp.exp(-A),yn*mp.exp(A),coeff[0]*e2+coeff[1]*e3)
                compare_forward(kind,source,Z,X,phi,64)
            # At t=0, one deck move exchanges the two nonnode axes.
            for a in map(mp.mpf,('-.2','.2')):
                x0,y0 = C.elementary_cusp_reference(0,a,nu)
                shifted = (mp.mpc(0),1/(nu*x0),-1) if a>0 else (
                    1/(nu*y0),mp.mpc(0),1)
                outputs = []
                for x,y,deck in ((x0,y0,0),shifted):
                    e2,e3,e4,a2,a3 = C._elementary_cusp_basis(0,inf)
                    c = C._elementary_real_coeffs(e4,e2,e3)
                    H = C.scalar_finite_height(x0,y0,nu)+deck
                    coeff = tuple(phi[j]+H*c[j] for j in range(2))
                    A = coeff[0]*a2+coeff[1]*a3
                    source = (x*mp.exp(-A),y*mp.exp(A),coeff[0]*e2+coeff[1]*e3)
                    outputs.append(C.elementary_global_forward(
                        kind,source,projection_steps=64,phase_steps=32,dps=55))
                error = max(abs(outputs[0][1]-outputs[1][1]),
                            *(abs(outputs[0][0][j]-outputs[1][0][j]) for j in range(2)),
                            *(abs(outputs[0][2][j]-outputs[1][2][j]) for j in range(2)))
                forward_error = max(forward_error,error)
                forward_cases += 1
                assert error<mp.mpf('1e-40')
        for X in map(mp.mpf,('.15','.5')):
            Z = tuple(mp.sqrt(1-X*X)*v for v in xi)
            kind,data = C.elementary_global_coordinates(
                Z,X,phi,projection_steps=64,phase_steps=32,dps=55)
            assert kind=='mid'
            z,w,t = data
            lt = mp.log(t)
            for dz,dw,dphi in ((2*mp.pi,0,(0,0)),
                               (-mp.j*t,1-t*t,(1,0)),
                               (-t,mp.j*(1+t*t),(0,1)),
                               (-mp.j*lt,2*t,(0,0))):
                shifted_phi = tuple(phi[j]+dphi[j] for j in range(2))
                compare_forward('mid',(z+dz,w+dw,t),Z,X,shifted_phi,64)
        # Independently chosen source points, compared after G(F(y)) in the
        # common middle chart modulo its four integer lattice directions.
        def as_middle(source_kind,source):
            if source_kind=='mid':
                return source
            x,y,w = source
            s = x*y
            if source_kind=='c0':
                t = s
                ell = mp.log(t)/(2*mp.j*mp.pi)
                return -mp.j*mp.log(y)+mp.pi*ell*(ell-1),w+2*t*ell,t
            t = 1/s
            wm = w*t*t
            return -mp.j*mp.log(y)+mp.j*wm/t,wm,t
        source_cases = []
        for inf,tc in ((False,mp.mpc('.02','.01')),
                       (True,mp.mpc('.08','-.03'))):
            y = mp.mpc('.4','.1')
            e2,e3,*_ = C._elementary_cusp_basis(tc,inf)
            source_cases.append(('ci' if inf else 'c0',
                                 (tc/y,y,mp.mpf('.13')*e2-mp.mpf('.09')*e3)))
        for t in (mp.mpc('.4','.15'),mp.mpc('.8','.3'),mp.mpc('2.3','-.4')):
            lt = mp.log(t)
            lattice = (mp.mpf('.12'),mp.mpf('.17'),mp.mpf('-.11'),mp.mpf('.34'))
            ez = (2*mp.pi,-mp.j*t,-t,-mp.j*lt)
            ew = (0,1-t*t,mp.j*(1+t*t),2*t)
            source_cases.append(('mid',
                                 (sum(lattice[j]*ez[j] for j in range(4)),
                                  sum(lattice[j]*ew[j] for j in range(4)),t)))
        for source_kind,source in source_cases:
            Z,X,angles = C.elementary_global_forward(
                source_kind,source,projection_steps=64,phase_steps=32,dps=55)
            back_kind,back = C.elementary_global_coordinates(
                Z,X,angles,projection_steps=64,phase_steps=32,dps=55)
            z,w,t = as_middle(source_kind,source)
            zb,wb,tb = as_middle(back_kind,back)
            lattice = C.scalar_lattice_coords(zb-z,wb-w,t,mp.log(t))
            error = max(abs(tb-t),*(abs(v-mp.nint(v)) for v in lattice))
            forward_error = max(forward_error,error)
            forward_cases += 1
            assert error<mp.mpf('1e-40')
        # The exact half-cell seam can return a different e4 deck
        # representative. Compare Y points in the common middle lattice.
        tc = mp.mpf('.1')*mp.exp(mp.j*mp.mpf('.7'))
        for inf in (False,True):
            nu = mp.exp(-2) if inf else mp.mpf(1)
            source_kind = 'ci' if inf else 'c0'
            for a in (mp.mpf('-.5'),mp.mpf('.5')):
                rho = C.elementary_cusp_logradius(a,tc,nu,dps=55)
                y = mp.exp(rho)
                source = (tc/y,y,mp.mpc(0))
                Z,X,angles = C.elementary_global_forward(
                    source_kind,source,projection_steps=512,phase_steps=32,dps=55)
                back_kind,back = C.elementary_global_coordinates(
                    Z,X,angles,projection_steps=512,phase_steps=32,dps=55)
                z,w,t = as_middle(source_kind,source)
                zb,wb,tb = as_middle(back_kind,back)
                lattice = C.scalar_lattice_coords(zb-z,wb-w,t,mp.log(t))
                error = max(abs(tb-t),*(abs(v-mp.nint(v)) for v in lattice))
                forward_error = max(forward_error,error)
                forward_cases += 1
                assert error<mp.mpf('1e-40')
        print(f'Elementary full forward: {forward_cases} reduced-count cap/free, '
              f'deck/lattice and independent source checks; error {float(forward_error):.1e}: PASS',
              flush=True)
        for inf in (False,True):
            nu = mp.exp(-2) if inf else mp.mpf(1)
            for a in map(mp.mpf,('-.499','-.2','.2','.499')):
                rad = C.inverse_height_moment(abs(a),0,0)
                x,y = (mp.sqrt(rad)/nu,mp.mpf(0)) if a>0 else (mp.mpf(0),mp.sqrt(rad))
                kind = 'ci' if inf else 'c0'
                e,norm = C.elliptic_bundle_embedding(kind,(x,y,a),dps=55,with_norm=True,elementary=True)
                assert mp.norm(mp.matrix(e))>mp.mpf('.999')
                _,h = C.elementary_restore_torus(kind,(x,y),phi)
                assert h[0]*h[1]==0 and (h[0]==0)==(x==0)
        print('Closed central fibres: 8 exact axis representatives, nonzero embedding and torus restoration: PASS',flush=True)

        # Projector identities and raw-base round trips including both central
        # circles and near-pole slab endpoints. Reconstruct u independently.
        raw_error = mp.mpf(0)
        projector_error = mp.mpf(0)
        projector_derivative = mp.mpf(0)
        rng = np.random.default_rng(20261008)
        for k in range(96):
            n = [mp.mpf(float(v)) for v in rng.normal(size=3)]
            n = [v/mp.sqrt(sum(w*w for w in n)) for v in n]
            X = list(map(mp.mpf,('-.999999999','-.999996','-.7','0','.6','.999996','.999999999','-.2')))[k%8]
            v = tuple((1-X*X)*w for w in n)
            kind,data = C.elliptic_raw_base_reference(v,X,dps=55,elementary=True)
            ref = C.elliptic_bundle_embedding(kind,data,dps=55,elementary=True)
            P = mp.matrix([[a*mp.conj(b) for b in ref] for a in ref])
            projector_error = max(projector_error,mp.norm(P*P-P),abs(sum(abs(a)**2 for a in ref)-1))
            if k<16:
                nn = mp.matrix(n)
                axis = mp.eye(3)[:,min(range(3),key=lambda j:abs(n[j]))]
                tangent = axis-(nn.T*axis)[0]*nn
                tangent /= mp.norm(tangent)
                second = mp.matrix([nn[1]*tangent[2]-nn[2]*tangent[1],
                                    nn[2]*tangent[0]-nn[0]*tangent[2],
                                    nn[0]*tangent[1]-nn[1]*tangent[0]])
                norms = []
                for direction in (tangent,second):
                    derivative = mp.diff(lambda h:C.elliptic_meridian_projector(
                        (nn+h*direction)/mp.norm(nn+h*direction),X,dps=55,elementary=True),0)
                    norms.append(mp.norm(derivative))
                derivative = mp.diff(lambda h:C.elliptic_meridian_projector(nn,X+h,dps=55,elementary=True),0)
                norms.append(mp.norm(derivative))
                projector_derivative = max(projector_derivative,mp.sqrt(sum(a*a for a in norms)))
            if kind=='mid':
                t,one,four,lt = data[:4]
                if len(data)==5:
                    inf = abs(t)>1
                    a = data[4]
                    tc = 1/t if inf else t
                else:
                    inf = abs(t)>1
                    tc = 1/t if inf else t
                    ss = lt/(2*mp.j*mp.pi)
                    xi4 = C.scalar_lattice_coords(mp.pi*ss*(ss-1),2*t*ss,t,lt)[3]
                    if inf:
                        lam = 1-C.scalar_flat_step((abs(t)-mp.mpf('1.2'))/mp.mpf('.4'))
                        a = -four-mp.mpf('.5')
                        aprime = -a+lam*(mp.arg(t)/(2*mp.pi)-xi4)+(1-lam)/2
                    else:
                        a = four-xi4-mp.mpf('.5')
            else:
                x,y,a = data
                inf = kind=='ci'
                tc = x*y
            A,B = (4*mp.pi)**2,mp.exp(4)
            if inf:
                ratio = tc*mp.sqrt((abs(tc)**2+A)/(A*abs(tc)**2+B))
                angle = aprime if kind=='mid' and len(data)==4 else -a+mp.mpf('.5')
                u1 = mp.exp(2*mp.pi*mp.j*angle)/mp.sqrt(1+abs(ratio)**2)
                u2 = ratio*u1
            else:
                ratio = tc*mp.sqrt((A+B*abs(tc)**2)/(1+A*abs(tc)**2))
                u2 = mp.exp(2*mp.pi*mp.j*a)/mp.sqrt(1+abs(ratio)**2)
                u1 = ratio*u2
            w = (u1.imag/(1+u1.real),u2.real/(1+u1.real)-1,u2.imag/(1+u1.real))
            rr = mp.sqrt(sum(a*a for a in w))
            backv = mp.matrix((2*w[0]/(1+rr*rr),-2*w[2]/(1+rr*rr),-2*w[1]/(1+rr*rr)))
            backX = (rr-1)/mp.sqrt(1+rr*rr)
            raw_error = max(raw_error,mp.norm(backv-mp.matrix(v)),abs(backX-X))
        assert raw_error<mp.mpf('1e-30') and projector_error<mp.mpf('1e-40')
        print(f'Closed raw base: 96 round trips/projectors; residual {float(raw_error):.1e}, projector {float(projector_error):.1e}: PASS',flush=True)
        assert projector_derivative<C.elementary_closed_parameters()['projector_bound']
        print(f'Projector derivative: 16 samples, measured Frobenius bound {float(projector_derivative):.6g}; proved operator bound 1e8000: PASS',flush=True)

        # Family logarithm: two chart expressions agree as REAL numbers,
        # despite many turns of phase and a moving interval endpoint.
        def gg(n,h):
            return mp.exp(mp.j*(12*n[0]+9*n[1]*n[2]+5*n[2]+40*h*n[0]+6*h*h))
        phase_error = mp.mpf(0)
        for k in range(12):
            theta = mp.mpf(k)*mp.pi/6
            nn = (mp.cos(theta),mp.sin(theta),mp.mpf(0))
            h = mp.mpf(k)/6-1
            expected = 12*nn[0]+40*h*nn[0]+6*h*h
            # The anchor north,h=0 has phase 5, whose chosen Arg is 5-2pi.
            expected -= 2*mp.pi
            a = C.elementary_phase_family(gg,nn,h,mp.mpf(0),1024,chart='north')
            b = C.elementary_phase_family(gg,nn,h,mp.mpf(0),1024,chart='south')
            phase_error = max(phase_error,abs(a-b),abs(a-expected))
        assert phase_error<mp.mpf('1e-40')
        print(f'Finite family log: 12 multi-turn real overlaps, error {float(phase_error):.1e}: PASS',flush=True)

        # The full pole/slab overlap, using reduced counts N=512,K=32.
        # Its proof for the actual formula uses the much larger fixed counts.
        xi = (mp.sqrt(mp.mpf('.7')),mp.sqrt(mp.mpf('.3'))*mp.exp(mp.j*mp.mpf('.7')))
        for sign in (-1,1):
            X = sign*(1-mp.mpf('4e-6'))
            Z = tuple(a*mp.sqrt(1-X*X) for a in xi)
            kind,cap = C.elementary_global_slice(Z,X,dps=55)
            diagnostic = {}
            other,slab = C.elementary_global_slice(Z,X,projection_steps=512,phase_steps=32,
                dps=55,force_slab=True,diagnostics=diagnostic)
            err = mp.norm(mp.matrix(cap)-mp.matrix(slab))
            assert kind==other and err<mp.mpf('1e-40')
            print(f'Finite N=512,K=32 cap overlap {sign:+d}: {float(err):.1e}; sampled denominators {dict((k,float(v)) for k,v in diagnostic.items())}: PASS',flush=True)

        # Actual local tensors from G, including exact poles and central axes.
        # First derivatives use arbitrary-precision differences, independently
        # of the generic SymPy naturality certificate C10.
        for sign in (-1,1):
            for vector in ((0,0,0,0),(mp.mpf('1e-5'),0,0,0),
                           (mp.mpf('1e-5'),mp.mpf('2e-5'),mp.mpf('-1e-5'),mp.mpf('1e-5'))):
                u = mp.matrix([*vector,*phi])
                J,DG = C.elementary_tensor_numeric(u,sign=sign,dps=55)
                assert abs(mp.det(DG))>mp.mpf('1e-6') and mp.norm(J*J+mp.eye(6))<mp.mpf('1e-40')
        # The replacement H_j differs from the historical blend only inside
        # this narrow shell. Exercise both directions and the actual cap tensor.
        shell_error = mp.mpf(0)
        shell_det = mp.inf
        for sign in (-1,1):
            inf = sign>0
            scale = mp.pi if inf else 2*mp.pi
            for radius in map(mp.mpf,('1.01e-8','1.25e-8','1.49e-8','2e-8')):
                local = (mp.sqrt(radius*mp.mpf('.7')),
                         mp.sqrt(radius*mp.mpf('.3'))*mp.exp(mp.j*mp.mpf('.7')))
                source = C.elementary_pole_patch(local,inf,dps=55)
                recovered = C.elementary_pole_forward(source,inf,dps=55)
                shell_error = max(shell_error,
                                  mp.norm(mp.matrix(recovered)-mp.matrix(local))/mp.sqrt(radius))
                if radius==mp.mpf('1.25e-8'):
                    h = scale*radius
                    X = sign*(1-h)/mp.sqrt(1+h*h)
                    U = (mp.matrix([[1,1],[1,-1]]) if inf else
                         mp.matrix([[1,1],[mp.j,-mp.j]]))/mp.sqrt(2)
                    Z = mp.sqrt(2*scale/(1+h*h))*U*mp.matrix(local)
                    assert abs(X*X+sum(abs(v)**2 for v in Z)-1)<mp.mpf('1e-50')
                    J,DG = C.elementary_tensor_numeric(
                        mp.matrix([Z[0].real,Z[0].imag,Z[1].real,Z[1].imag,*phi]),
                        sign=sign,dps=55)
                    shell_det = min(shell_det,abs(mp.det(DG)))
                    assert mp.norm(J*J+mp.eye(6))<mp.mpf('1e-40')
        assert shell_error<mp.mpf('1e-40') and shell_det>mp.mpf('1e-6')
        print(f'Elementary pole shells: 8 forward/inverse cases, relative error {float(shell_error):.1e}; '
              f'2 cap tensors, minimum |det DG|={float(shell_det):.6g}: PASS',flush=True)
        # The radial conjugation is flatter than the precision used by the
        # raw stereographic formulas very close to its inner and outer seams.
        seam_error = mp.mpf(0)
        for inf in (False,True):
            for radius in [mp.mpf('1e-8')*(1+gap) for gap in
                           map(mp.mpf,('1e-40','1e-10','1e-6','5.5e-5','6e-5','7e-5','1e-4'))]+[
                           mp.mpf('1.5e-8')*(1-mp.mpf('1e-6')),
                           mp.mpf('1.5e-8')*(1-mp.mpf('1e-30'))]:
                local = (mp.sqrt(radius*mp.mpf('.7')),
                         mp.sqrt(radius*mp.mpf('.3'))*mp.exp(mp.j*mp.mpf('.7')))
                source = C.elementary_pole_patch(local,inf,dps=55)
                recovered = C.elementary_pole_forward(source,inf,dps=55)
                seam_error = max(seam_error,
                                 mp.norm(mp.matrix(recovered)-mp.matrix(local))/mp.sqrt(radius))
        assert seam_error<mp.mpf('1e-40')
        print(f'Elementary pole seams: 18 near-inner/outer forward/inverse cases, '
              f'relative error {float(seam_error):.1e}: PASS',flush=True)
        # A free point exercises a full projection product, with no cap cutoff.
        X = mp.mpf('.4')
        Z = tuple(a*mp.sqrt(1-X*X) for a in xi)
        u = mp.matrix([Z[0].real,Z[0].imag,Z[1].real,Z[1].imag,*phi])
        J,DG = C.elementary_tensor_numeric(u,projection_steps=512,phase_steps=32,dps=55)
        assert abs(mp.det(DG))>mp.mpf('1e-6') and mp.norm(J*J+mp.eye(6))<mp.mpf('1e-40')
        q = [Z[0].real,Z[0].imag,Z[1].real,Z[1].imag,X]
        alternate = mp.matrix([*q[1:],*phi])
        J_alt,_ = C.elementary_tensor_numeric(alternate,missing=0,projection_steps=512,phase_steps=32,dps=55)
        transition = mp.matrix(6)
        for i in range(3):
            transition[i,i+1] = 1
        for j in range(4):
            transition[3,j] = -q[j]/X
        transition[4,4] = transition[5,5] = 1
        assert mp.norm(J_alt-transition*J*transition**-1)<mp.mpf('1e-35')
        shifted = u+mp.matrix([0,0,0,0,1,-1])
        J_periodic,_ = C.elementary_tensor_numeric(shifted,projection_steps=512,phase_steps=32,dps=55)
        assert mp.norm(J_periodic-J)<mp.mpf('1e-35')
        print(f'Local tensors: 6 pole/axis points and one free N=512 point; free |det DG|={float(abs(mp.det(DG))):.6g}: PASS',flush=True)
        print('New free tensor: round-sphere chart transition and integer torus shifts agree: PASS',flush=True)
        for height in (mp.mpf(0),1/mp.sqrt(5)):
            # Hopf direction north; these are the two a=1/2 central-circle
            # points u=(1,0) and u=(0,-1) in the orbit S3, respectively.
            central = mp.matrix([mp.sqrt(1-height*height)*mp.sin(mp.mpf('.37')),0,0,height,*phi])
            Jc,Mc = C.elementary_tensor_numeric(central,missing=0,projection_steps=512,phase_steps=32,dps=55)
            assert abs(mp.det(Mc))>mp.mpf('1e-6') and mp.norm(Jc*Jc+mp.eye(6))<mp.mpf('1e-40')
            print(f'New central-circle tensor X={float(height):.6f}: |det DG|={float(abs(mp.det(Mc))):.6g}, transverse derivatives included: PASS',flush=True)
    print('ALL PASS (closed blocks; reduced-count checks are not a global numerical certificate)',flush=True)


def check_elementary_tensor_differences():
    """Independent finite-difference N and a nonintegrable control for new G."""
    import mpmath as mp
    from extra.s4t2_explicit import cusp as C
    with mp.workdps(40):
        basis = mp.eye(6)
        def nij(J,derivatives):
            terms = [[],[],[]]
            for k in range(6):
                for i in range(6):
                    for j in range(6):
                        terms[0].append(sum(J[l,i]*derivatives[l][k,j] for l in range(6)))
                        terms[1].append(-sum(J[l,j]*derivatives[l][k,i] for l in range(6)))
                        terms[2].append(-sum(J[k,l]*(derivatives[i][l,j]-derivatives[j][l,i]) for l in range(6)))
            vectors = [mp.matrix(v) for v in terms]
            return mp.norm(sum(vectors,mp.matrix(216,1)))/max(1,sum(mp.norm(v) for v in vectors))
        def perturb(J,u):
            E = mp.matrix(6)
            # Mix a sphere direction with a torus direction. A shear inside
            # one complex line can remain integrable and is not a valid control.
            E[0,4] = mp.mpf('.2')*mp.sin(u[0])
            return (basis+E)*J*(basis-E)
        freeZ = (mp.sqrt(mp.mpf('.7')*mp.mpf('.84')),
                 mp.sqrt(mp.mpf('.3')*mp.mpf('.84'))*mp.exp(mp.j*mp.mpf('.7')))
        points = [(-1,mp.matrix([mp.mpf('1e-5'),mp.mpf('2e-5'),mp.mpf('-1e-5'),mp.mpf('1e-5'),mp.mpf('.17'),mp.mpf('-.23')])),
                  (1,mp.matrix([mp.mpf('1e-5'),0,0,mp.mpf('1e-5'),mp.mpf('.17'),mp.mpf('-.23')])),
                  (1,mp.matrix([freeZ[0].real,freeZ[0].imag,freeZ[1].real,freeZ[1].imag,mp.mpf('.17'),mp.mpf('-.23')]))]
        for index,(sign,u) in enumerate(points):
            def tensor(v):
                return C.elementary_tensor_numeric(v,sign=sign,projection_steps=512,phase_steps=32,dps=40)[0]
            J = tensor(u)
            results = []
            for h in (mp.mpf('2e-7'),mp.mpf('1e-7')):
                derivatives,controls = [],[]
                for j in range(6):
                    up,um = u+h*basis[:,j],u-h*basis[:,j]
                    plus,minus = tensor(up),tensor(um)
                    derivatives.append((plus-minus)/(2*h))
                    controls.append((perturb(plus,up)-perturb(minus,um))/(2*h))
                value,control = nij(J,derivatives),nij(perturb(J,u),controls)
                assert value<mp.mpf('1e-5') and control>mp.mpf('1e-5'), (index,float(h),float(value),float(control))
                results.append(float(value))
                print(f'New tensor FD N point {index}, h={float(h):.1e}: normalized N={float(value):.2e}, nonintegrable control={float(control):.2e}: PASS',flush=True)
            # A central-difference truncation should fall quadratically; allow
            # roundoff domination only far below the assertion threshold.
            assert results[1]<.3*results[0] or results[1]<1e-25
    print('ALL PASS (new finite-instance tensor FD checks)',flush=True)


if __name__ == '__main__':
    if '--closed-nijenhuis' in sys.argv:
        check_elementary_tensor_differences()
        raise SystemExit(0)
    if '--closed' in sys.argv:
        check_elementary_closed_blocks()
        raise SystemExit(0)
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
    torch.manual_seed(5)
    q = NR.S4T2.sample(n)
    t0 = time.time()
    tt, lg, x = TI.initial_guess(q)
    guess_err = 0.0
    for key, (mask, P) in TI.to_chart(tt, lg, x).items():
        if mask.any():
            guess_err = max(guess_err, (TI.CHART_F[key](P[mask]) - q[mask]).abs().max().item())
    sol = TI.Finv(q)
    res = torch.cat([v[2] for v in sol.values()])
    print(f'{n} uniform points of S^4 x T^2: charts used ' + ', '.join(f'{k}: {len(v[0])}' for k, v in sol.items()))
    print(f'   explicit chain (before Newton): max |F(p0) - q| = {guess_err:.1e}')
    print(f'   after Newton polish: max |F(p) - q| = {res.max():.1e}, median {res.median():.1e}  ({time.time() - t0:.0f} s)')
    worst = []
    for key, (idx, p, r) in sol.items():
        qq, A, H = NR.DF_and_D2F(TI.CHART_F[key], p)
        J, DJ = NR.transport(A, H, qq)
        Nn, scale, Pq = NR.nij(qq, J, DJ)
        worst.append((Nn / scale).max().item())
        print(f'   {key}: J at F^-1(q): |J^2 + P| max {(J @ J + Pq).flatten(1).norm(dim=-1).max():.1e}; '
              f'|N| / scale max {(Nn / scale).max():.1e}; |J| median {J.flatten(1).norm(dim=-1).median():.1f}')
    ok = res.max() < 1e-10
    print('PASS' if ok else 'FAIL')
    assert ok

    # Exact poles used to return no inverse entries at all, silently producing J=0.
    # Check chart coverage as well as residuals; include points on each side of the
    # new closed-inverse boundary, and a mixed batch exercising merged chart entries.
    from extra.s4t2_explicit import tlift as L, tfull as Fm
    gen = torch.Generator().manual_seed(73)
    targets = []
    for sign, scale in ((-1, L.S0), (1, L.SINF)):
        radii = torch.tensor([0., 1e-20, L.EPS1 / 2, L.EPS1 * (1 - 1e-6), L.EPS1 * (1 + 1e-6), L.EPS2])
        radii = radii.repeat_interleave(8)
        Z = torch.randn(len(radii), 4, generator=gen)
        Z = Z / Z.norm(dim=-1, keepdim=True) * torch.sqrt(scale * radii[:, None])
        height = sign * torch.sqrt(1 - scale * radii)
        phi = torch.rand(len(radii), 2, generator=gen) - 0.5
        targets.append(torch.cat([Z, height[:, None], Fm.torus_embed(phi)], -1))
    targets = torch.cat(targets)
    closed_sol = TI.Finv(targets)
    indices = torch.cat([v[0] for v in closed_sol.values()])
    assert torch.equal(indices.sort().values, torch.arange(len(targets)))
    ball_error = 0.0
    for key, (idx, p, residual) in closed_sol.items():
        assert torch.isfinite(p).all() and torch.isfinite(residual).all()
        error = (TI.CHART_F[key](p) - targets[idx]).abs().max().item()
        assert error < 1e-10
        ball_error = max(ball_error, error)
        # The closed piece itself needs no Newton, including at the exact pole.
        scale = L.S0 if key == 'c0' else L.SINF
        inner = (targets[idx, :4] ** 2).sum(-1) < scale * L.EPS1
        assert torch.allclose(p[inner], TI.linear_ball_inverse(targets[idx][inner], key), atol=1e-14, rtol=1e-14)
    print(f'Linear balls: {len(targets)} pole/near-pole/boundary points, complete chart coverage, max residual {ball_error:.1e}: PASS')
    check_transport_base()
    check_transport_base(algebraic_radial=True)
    check_transport_base(algebraic_radial=True,algebraic_profile=True)
    check_algebraic_axis_jets()
