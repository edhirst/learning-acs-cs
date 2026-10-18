"""Focused executable verification of direct.py.

Exact identities are distinguished from numerical comparisons. The written
domain/smoothness contracts remain mathematical inputs for human review.
Run from the repository root:
  .venv/bin/python3 -B extra/s4t2_explicit/sym/direct_certificate.py
"""
from pathlib import Path
import ast
import inspect
import sys
import time
import types
import mpmath as mp
import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from extra.s4t2_explicit import direct as D


def check(name, condition):
    if not condition:
        raise AssertionError(name)
    print(name + ": PASS", flush=True)


def exact_frame():
    """All five sphere-coordinate choices, with the constraint ||x||^2=1."""
    x = sp.Matrix(sp.symbols('x0:5', real=True))
    c = D.frame_brackets(list(x),dtype=object)
    for i in range(5):
        for j in range(5):
            Vi,Vj = sp.eye(5)[:,i]-x[i]*x,sp.eye(5)[:,j]-x[j]*x
            bracket = Vj.jacobian(x)*Vi-Vi.jacobian(x)*Vj
            assert all(sp.expand(v)==0 for v in bracket-x[i]*Vj+x[j]*Vi)
            assert sp.expand(Vi[j]-Vj[i])==0
            assert all(c[k,i,j]==x[i]*int(k==j)-x[j]*int(k==i) for k in range(6))
    assert all(v==0 for v in c[:,:,5].flat) and all(v==0 for v in c[:,5,:].flat)
    ideal = sp.groebner([x.dot(x)-1], *x)
    for missing in range(5):
        keep = [i for i in range(5) if i != missing]
        T, inverse = map(sp.Matrix,D.frame_matrices(list(x),missing,dtype=object))
        derivatives = D.frame_derivatives(list(x),missing,dtype=object)
        for k,i in enumerate(keep):
            expected = T.diff(x[i])-x[i]/x[missing]*T.diff(x[missing])
            assert all(sp.cancel(v)==0 for v in expected-sp.Matrix(derivatives[k]))
        assert not any(derivatives[4:].flat)
        for residual in (T*inverse-sp.eye(6), inverse*T-sp.eye(6)):
            for value in residual:
                numerator = sp.fraction(sp.cancel(value))[0]
                assert ideal.reduce(numerator)[1] == 0
    check("D1 exact global frame, brackets, derivatives and all hemisphere inverses", True)


def exact_tensor():
    """Check A*N without normalizing A or expanding a symbolic 6x6 inverse.

    From A J=J0 A and A partial_l J=J0 H_l-H_l J, multiplying
    the coordinate N formula by A gives the polynomial below. Each
    coefficient vanishes for every J and every symmetric Hessian H.
    Invertibility of A is supplied by the explicit map/inverse contract.
    """
    J0 = sp.diag(*([sp.Matrix([[0, -1], [1, 0]])]*3))
    J = sp.Matrix(6, 6, sp.symbols('j0:36'))
    H = [[[sp.Symbol(f'h{a}_{min(i,j)}_{max(i,j)}')
           for j in range(6)] for i in range(6)] for a in range(6)]
    # Bind the actual symbolic tensor constructor to the formulas used below.
    # Only its derivative inputs are formal; its implementation is executed.
    A = sp.ImmutableMatrix(6,6,sp.symbols('a0:36'))
    class FormalTensor(D.SymbolicTensor):
        def differentiate(self):
            return A,tuple(sp.Matrix(h) for h in H),()
    tensor = FormalTensor(types.SimpleNamespace(sp=sp),(),'mid',4,1)
    actual_J,actual_dJ,_ = tensor.tensor_expressions()
    inverse = sp.Inverse(A)
    expected_J = sp.MatMul(inverse,J0,A,evaluate=False)
    assert actual_J==expected_J
    for k in range(6):
        Hk = sp.ImmutableMatrix(6,6,lambda a,i:H[a][i][k])
        expected_dJ = sp.MatMul(inverse,sp.MatAdd(
            sp.MatMul(J0,Hk,evaluate=False),
            sp.MatMul(-1,Hk,expected_J,evaluate=False),evaluate=False),evaluate=False)
        assert actual_dJ[k]==expected_dJ
    # Exact square calculation in the associative word algebra, A B=B A=I.
    def reduce_word(word):
        word = list(word)
        changed = True
        coefficient = 1
        while changed:
            changed = False
            for i in range(len(word)-1):
                pair = word[i:i+2]
                if pair in (['A','B'], ['B','A'], ['J0','J0']):
                    if pair == ['J0','J0']:
                        coefficient *= -1
                    del word[i:i+2]
                    changed = True
                    break
        return coefficient, word
    assert J0*J0 == -sp.eye(6)
    assert reduce_word(['B','J0','A','B','J0','A']) == (-1, [])
    count = 0
    for a in range(6):
        for i in range(6):
            for j in range(6):
                expression = H[a][j][i]-H[a][i][j]
                expression += sum(J0[a,b]*(
                    J[l,i]*H[b][j][l]-J[l,j]*H[b][i][l]
                    +H[b][l][i]*J[l,j]-H[b][l][j]*J[l,i])
                    for b in range(6) for l in range(6))
                expression += sum(-J[l,i]*H[a][m][l]*J[m,j]
                                  +J[l,j]*H[a][m][l]*J[m,i]
                                  for l in range(6) for m in range(6))
                assert sp.expand(expression) == 0
                count += 1
    # Independent nonintegrable control: symmetric Hessians are essential.
    dJ = [sp.zeros(6) for _ in range(6)]
    K = sp.zeros(6); K[0, 0] = 1
    dJ[2] = J0*K-K*J0
    bad = D.nijenhuis_components(np.array(J0,dtype=object),np.array(dJ,dtype=object))
    assert any(sp.expand(v)!=0 for v in bad.flat)
    good_dJ = []
    for k in range(6):
        Hk = sp.Matrix(6,6,lambda a,i:H[a][i][k])
        good_dJ.append(J0*Hk-Hk*J0)
    good = D.nijenhuis_components(np.array(J0,dtype=object),np.array(good_dJ,dtype=object))
    assert all(sp.expand(v)==0 for v in good.flat)
    check(f"D2 exact square and {count} Nijenhuis component reductions; negative control", True)


def exact_phase():
    u, v = sp.symbols('u v', real=True)
    du, dv, eu, ev, hu, hv = sp.symbols('du dv eu ev hu hv', real=True)
    a, b = sp.symbols('a b', real=True)
    U = u+du*a+eu*b+hu*a*b
    V = v+dv*a+ev*b+hv*a*b
    g, ga, gb, gab = u+sp.I*v, du+sp.I*dv, eu+sp.I*ev, hu+sp.I*hv
    angle = sp.atan2(V, U)
    at = {a: 0, b: 0}
    assert sp.simplify(sp.diff(angle,a).subs(at)-sp.im(ga/g)) == 0
    assert sp.simplify(sp.diff(angle,a,b).subs(at)-sp.im(gab/g-ga*gb/g**2)) == 0
    # Quotient logarithmic derivatives telescope in the actual finite sum.
    x, y, dx, dy = sp.symbols('x y dx dy', nonzero=True)
    quotient = (x+a*dx)/(y+a*dy)
    assert sp.cancel((sp.diff(quotient,a)/quotient).subs(a,0)-dx/x+dy/y) == 0
    assert sp.Rational(22,7)*(1900+2700)/D.PHASE_STEPS < sp.Rational(1,4)
    assert 1-sp.Rational(1,4)**2/2 == sp.Rational(31,32)
    check("D3 exact finite-phase derivative rule and increment margin (given the written derivative bounds)", True)


def exact_jet_rules(unary_override=None):
    """Execute the actual numeric unary method bodies with symbolic scalars."""
    z = sp.Symbol('z', positive=True)
    class Capture(D.Jet):
        def __init__(self, value):
            self.value = value
        def unary(self, value, first, second):
            return value, first, second
    symbolic_mp = types.SimpleNamespace(
        exp=sp.exp, expm1=lambda x:sp.exp(x)-1, log=sp.log,
        log1p=lambda x:sp.log(1+x), sin=sp.sin, cos=sp.cos,
        sinh=sp.sinh, cosh=sp.cosh, asinh=sp.asinh, sqrt=sp.sqrt)
    for name in ('exp','expm1','log','log1p','sin','cos','sinh','cosh','asinh'):
        original = (unary_override or {}).get(name,getattr(D.Arithmetic, name))
        environment = dict(original.__globals__, mp=symbolic_mp)
        actual = types.FunctionType(original.__code__, environment)
        value, first, second = actual(Capture(z))
        assert sp.simplify(sp.diff(value,z)-first) == 0
        assert sp.simplify(sp.diff(first,z)-second) == 0
    # Symmetric composition Hessian, using independent first derivatives.
    fu, fv, fuu, fuv, fvv = sp.symbols('fu fv fuu fuv fvv')
    ui, uj, vi, vj, uij, vij = sp.symbols('ui uj vi vj uij vij')
    mixed = fu*uij+fv*vij+fuu*ui*uj+fuv*(ui*vj+vi*uj)+fvv*vi*vj
    assert sp.expand(mixed-mixed.xreplace({ui:uj,uj:ui,vi:vj,vj:vi})) == 0
    # Bind the actual product and unary rule; no digest is called a proof.
    def returned(method):
        import textwrap
        tree = ast.parse(textwrap.dedent(inspect.getsource(method)))
        return ast.dump([n for n in tree.body[0].body if isinstance(n,ast.Return)][-1].value)
    def expression(text):
        return ast.dump(ast.parse(text, mode='eval').body)
    assert returned(D.Jet.unary) == expression(
        'Jet(value, first*self.gradient, first*self.hessian+second*np.outer(self.gradient,self.gradient))')
    assert returned(D.Jet.__mul__) == expression(
        'Jet(self.value*other.value, self.gradient*other.value+self.value*other.gradient, '
        'self.hessian*other.value+self.value*other.hessian+'
        'np.outer(self.gradient,other.gradient)+np.outer(other.gradient,self.gradient))')
    assert returned(D.Jet.__pow__) == expression(
        'self.unary(self.value**power, power*self.value**(power-1), '
        'power*(power-1)*self.value**(power-2))')
    p = sp.Symbol('p',real=True)
    assert sp.simplify(sp.diff(z**p,z,2)-p*(p-1)*z**(p-2))==0
    assert returned(D.phase_jet) == expression(
        'Jet(value, imag_array(g.gradient/g.value), '
        'imag_array(g.hessian/g.value-np.outer(g.gradient,g.gradient)/g.value**2))')
    assert returned(D.Arithmetic.arg) == expression(
        'Jet(mp.arg(x.value), imag_array(x.gradient/x.value), '
        'imag_array(x.hessian/x.value-np.outer(x.gradient,x.gradient)/x.value**2))')
    check("D4 actual unary derivative implementations and finite product/chain rules", True)


def symbolic_primitive_rules():
    x = sp.Symbol('x',positive=True)
    y = sp.Symbol('y',real=True)
    def unpack(number):
        substitution = {}
        for node in number.graph.assignments:
            substitution[node.symbol] = node.expression.xreplace(substitution)
        return (number.re+sp.I*number.im).xreplace(substitution)
    for name in ('exp','expm1','log','log1p','sin','cos','arg','sinh','cosh','asinh'):
        graph = D.ExpressionGraph()
        real_only = name in ('sinh','cosh','asinh')
        z = graph.number(mp.mpc(2,0 if real_only else .3),x,0 if real_only else y)
        actual = unpack(D.symbolic_unary(name,z))
        arg = x if real_only else x+sp.I*y
        if name=='expm1':
            expected = sp.exp(arg)-1
        elif name=='log1p':
            expected = sp.log(1+arg)
        elif name=='arg':
            expected = sp.atan2(y,x)
        else:
            expected = getattr(sp,name)(arg)
        expected = sp.expand_complex(expected).replace(
            lambda e:e.func==sp.arg,
            lambda e:sp.atan2(sp.im(e.args[0]),sp.re(e.args[0])))
        assert sp.simplify(sp.expand_complex(actual-expected))==0,name
    for operation,expected in [
            (lambda z:z+z,2*(x+sp.I*y)),
            (lambda z:z*z,(x+sp.I*y)**2),
            (lambda z:1/z,1/(x+sp.I*y)),
            (lambda z:z.conjugate(),x-sp.I*y),
            (lambda z:D.norm_squared(z),x*x+y*y)]:
        graph = D.ExpressionGraph()
        z = graph.number(mp.mpc(2,.3),x,y)
        assert sp.simplify(sp.expand_complex(unpack(operation(z))-expected))==0
    # Check the actual named derivative compiler against direct SymPy
    # differentiation of one fully expanded, mixed six-variable expression.
    graph = D.ExpressionGraph()
    q = [graph.number(mp.mpf(i+1)/10,s,0) for i,s in enumerate(graph.inputs)]
    u = q[0]*q[1]+q[2]
    v = q[3]-q[4]*q[5]
    z = D.symbolic_unary('sin',u)*D.symbolic_unary('exp',v)
    z = z+D.symbolic_unary('log1p',u*u+v*v)
    A,H,program = graph.differentiate((z.re,))
    definitions = {}
    for node in graph.active_assignments((z.re,)):
        definitions[node.symbol] = node.expression.xreplace(definitions)
    exact = z.re.xreplace(definitions)
    for symbol,expression in program:
        definitions[symbol] = expression.xreplace(definitions)
    for i,qi in enumerate(graph.inputs):
        assert sp.simplify(A[0,i].xreplace(definitions)-sp.diff(exact,qi))==0
        for j,qj in enumerate(graph.inputs):
            assert sp.simplify(H[0][i,j].xreplace(definitions)-sp.diff(exact,qi,qj))==0
    # Mutation control: a wrong second derivative in the actual unary kernel
    # must be rejected, without altering the imported production module.
    import textwrap
    source = textwrap.dedent(inspect.getsource(D.Arithmetic.sin))
    source = source.replace('@staticmethod\n','').replace('-mp.sin(x.value)','mp.sin(x.value)')
    namespace = dict(vars(D))
    exec(compile(source,'<negative derivative control>','exec'),namespace)
    try:
        exact_jet_rules({'sin':namespace['sin']})
    except AssertionError:
        pass
    else:
        raise AssertionError('incorrect second derivative escaped the verifier')
    check("D10 exact symbolic arithmetic agrees with elementary functions; derivative mutation rejected",True)


def formula_correspondence():
    """Compare complete executable bodies, with an explicit small rewrite.

    The original checks are not executed or imported. Their formulas and
    cusp.py are read as source. This test establishes correspondence to the
    inherited proofs, not the truth of those written proofs.
    """
    import copy
    root = Path(D.__file__).parent
    old_tree = ast.parse((root/'sym/chart_identities.py').read_text())
    cusp_tree = ast.parse((root/'cusp.py').read_text())
    old = {n.name:n for n in old_tree.body if isinstance(n,ast.FunctionDef)}
    cusp = {n.name:n for n in cusp_tree.body if isinstance(n,ast.FunctionDef)}
    names = {
        'direct_middle_matrix58':'middle_lift','direct_cusp_matrix59':'cusp_lift',
        'collar_offset61':'collar_offset','extended_cusp_matrix62':'extended_cusp_lift',
        'boundary_exact66':'boundary_mismatch','direct_global_section65':'circle_section',
        'direct_global_coordinates65':'coordinates','direct_global_forward65':'forward',
        'elliptic_raw_base_reference':'base_reference',
        'elliptic_bundle_embedding':'bundle_embedding',
        'elliptic_pole_raw_inverse':'pole_raw_inverse','elementary_global_slice':'pole_slice',
    }
    reverse = {v:k for k,v in names.items()}
    methods = {n.name:n for n in ast.parse(inspect.getsource(D.Formulas)).body[0].body
               if isinstance(n,ast.FunctionDef) and n.name != '__init__'}
    source_names = {reverse.get(n,n) for n in methods}

    class Rewrite(ast.NodeTransformer):
        def __init__(self, local_names):
            self.local_names = local_names
        def visit_Name(self,n):
            if n.id == 'elementary':
                return ast.copy_location(ast.Constant(True),n)
            if n.id in source_names and n.id not in self.local_names:
                return ast.Attribute(ast.Name('self',ast.Load()),names.get(n.id,n.id),ast.Load())
            if n.id == 'boundary_phase64':
                return ast.Name('phase_value',ast.Load())
            return n
        def visit_Attribute(self,n):
            if isinstance(n.value,ast.Name) and n.value.id == 'cusp30':
                return ast.Attribute(ast.Name('self',ast.Load()),names.get(n.attr,n.attr),ast.Load())
            return self.generic_visit(n)
        def visit_Import(self,n):
            if len(n.names)==1 and n.names[0].name=='mpmath':
                return None
            return n
        def visit_Call(self,n):
            n = self.generic_visit(n)
            n.keywords = [k for k in n.keywords if k.arg!='elementary']
            return n
        def visit_FunctionDef(self,n):
            n = self.generic_visit(n)
            if n.args.args and n.args.args[-1].arg=='elementary':
                n.args.args.pop()
                n.args.defaults.pop()
            return n
        def visit_If(self,n):
            if ast.unparse(n.test)=='_elementary_precision_shortcuts.get()':
                return None
            n = self.generic_visit(n)
            if isinstance(n.test,ast.Constant):
                return n.body if n.test.value else n.orelse
            return n
        def visit_IfExp(self,n):
            n = self.generic_visit(n)
            if isinstance(n.test,ast.Constant):
                return n.body if n.test.value else n.orelse
            return n
        def visit_BoolOp(self,n):
            n = self.generic_visit(n)
            if isinstance(n.op,ast.And):
                if any(isinstance(v,ast.Constant) and not v.value for v in n.values):
                    return ast.Constant(False)
                n.values = [v for v in n.values if not (isinstance(v,ast.Constant) and v.value)]
            return n.values[0] if len(n.values)==1 else n
        def visit_BinOp(self,n):
            n = self.generic_visit(n)
            if (isinstance(n.op,ast.Pow) and isinstance(n.left,ast.Call)
                and isinstance(n.left.func,ast.Name) and n.left.func.id=='abs'
                and isinstance(n.right,ast.Constant) and isinstance(n.right.value,int)
                and n.right.value>0 and n.right.value%2==0):
                square = ast.Call(ast.Name('norm_squared',ast.Load()),n.left.args,[])
                return square if n.right.value==2 else ast.BinOp(square,ast.Pow(),ast.Constant(n.right.value//2))
            return n

    def without_docstrings(node):
        class Strip(ast.NodeTransformer):
            def visit_FunctionDef(self,n):
                n = self.generic_visit(n)
                if n.body and isinstance(n.body[0],ast.Expr) and isinstance(n.body[0].value,ast.Constant):
                    n.body = n.body[1:]
                return n
        return ast.dump(Strip().visit(copy.deepcopy(node)),include_attributes=False)

    for name, actual in methods.items():
        original_name = reverse.get(name,name)
        original = copy.deepcopy((old if original_name in old else cusp)[original_name])
        if name == 'pole_slice':
            with_node = next(n for n in original.body if isinstance(n,ast.With))
            index = next(i for i,n in enumerate(with_node.body) if isinstance(n,ast.If)
                         and isinstance(n.test,ast.BoolOp))
            with_node.body = with_node.body[:index]+with_node.body[index].body
            original.body = [with_node]
        local_names = {n.id for n in ast.walk(original)
                       if isinstance(n,ast.Name) and isinstance(n.ctx,ast.Store)}
        local_names |= {n.arg for n in ast.walk(original) if isinstance(n,ast.arg)}
        local_names |= {n.name for n in ast.walk(original) if isinstance(n,ast.FunctionDef) and n is not original}
        expected = Rewrite(local_names).visit(original)
        if name=='pole_slice':
            with_node = expected.body[0]
            with_node.body = [n for n in with_node.body if not
                             (isinstance(n,ast.Assign) and isinstance(n.targets[0],ast.Name)
                              and n.targets[0].id in ('alpha','distance'))]
            expected.args = ast.parse('def f(Z,X,*,dps=50): pass').body[0].args
        if name == 'base_reference':
            with_node = next(n for n in expected.body if isinstance(n,ast.With))
            outer = next(n for n in with_node.body if isinstance(n,ast.If) and isinstance(n.test,ast.BoolOp))
            outer.body = [n for n in outer.body if not (isinstance(n,ast.If) and ast.unparse(n.test)=='tc == 0')]
        expected.name = name
        expected.args.args.insert(0,ast.arg('self'))
        index = 1 if expected.body and isinstance(expected.body[0],ast.Expr) else 0
        expected.body.insert(index,ast.parse('mp = self.ops').body[0])
        assert without_docstrings(expected) == without_docstrings(actual), name
    contracts = D.proof_contracts()
    covered = [name for contract in contracts for name in contract.methods]
    assert len(covered) == len(set(covered)) and set(covered)==set(methods)
    # Vectorized phase is an approximate evaluation of the same fixed sum.
    for old_name, new_name in [('boundary_transition62','boundary_transition'),('boundary_phase64','phase_value')]:
        expected = copy.deepcopy(old[old_name])
        class Rename(ast.NodeTransformer):
            def visit_Name(self,n):
                n.id = {'np57':'np','boundary_transition62':'boundary_transition'}.get(n.id,n.id)
                return n
        expected = Rename().visit(expected)
        expected.body = [n for n in expected.body if not
                         (isinstance(n,ast.Assign) and ast.unparse(n)=='np = np')]
        expected.name = new_name
        actual = ast.parse(inspect.getsource(getattr(D,new_name))).body[0]
        assert without_docstrings(expected)==without_docstrings(actual), new_name
    assert D.PHASE_STEPS==inspect.signature(D.phase_value).parameters['steps'].default==65536
    check(f"D6 complete correspondence of {len(methods)} formula bodies; {len(contracts)} written contracts", True)


def round_point(u1,u2):
    p = mp.matrix([u1.imag,-u2.imag,1+u1.real-u2.real])
    n = p/mp.norm(p)
    beta = mp.norm(p)/(1+u1.real)
    X = (beta-1)/mp.sqrt(1+beta*beta)
    xi = D.Formulas()._elementary_hopf_reference(n)
    Z = [mp.sqrt(1-X*X)*v*mp.exp(mp.j*mp.mpf('.37')) for v in xi]
    return [Z[0].real,Z[0].imag,Z[1].real,Z[1].imag,X,mp.mpf('.17'),mp.mpf('-.23')]


def radial_point(radius,theta=mp.mpf('.81'),phi=mp.mpf('-.63')):
    ratio = radius**2*((4*mp.pi)**2+mp.exp(4)*radius**2)/(1+(4*mp.pi)**2*radius**2)
    a,b = mp.sqrt(ratio/(1+ratio)),1/mp.sqrt(1+ratio)
    return round_point(a*mp.exp(mp.j*(theta+phi)),b*mp.exp(mp.j*phi))


def chart_input(point,missing=None):
    if missing is None:
        missing = max(range(5),key=lambda i:abs(point[i]))
    return [point[i] for i in range(5) if i!=missing]+point[5:],missing,1 if point[missing]>0 else -1


def evaluate_symbolic_jets(symbolic,A,H,assignments,values):
    """Numerical observation using the independently differentiated exact DAG."""
    from sympy.printing.pycode import MpmathPrinter
    class Printer(MpmathPrinter):
        def _print_Function(self, expr):
            if expr.func.__name__=='FinitePhase65536':
                return 'finite_phase('+','.join(self._print(a) for a in expr.args)+')'
            return super()._print_Function(expr)
    printer = Printer()
    lines = ['def evaluate('+','.join(map(str,symbolic.graph.inputs))+'):']
    for node in symbolic.graph.active_assignments(symbolic.outputs):
        lines.append('    '+str(node.symbol)+' = '+printer.doprint(node.expression))
    for symbol,expression in assignments:
        lines.append('    '+str(symbol)+' = '+printer.doprint(expression))
    flat = list(A)+[v for h in H for v in h]
    lines.append('    return ('+','.join(printer.doprint(e) for e in flat)+',)')
    env = {'mpmath':mp,'finite_phase':lambda th,ph:mp.mpf(D.phase_value(float(th),float(ph))[0])}
    exec(compile('\n'.join(lines),'<exact derivative program>','exec'),env)
    result = env['evaluate'](*values)
    return np.array(result[:36],float).reshape(6,6),np.array(result[36:],float).reshape(6,6,6)


def symbolic_cases():
    fixtures = [
        ('south pole', [0,0,0,0,sp.Rational(17,100),-sp.Rational(23,100)],4,-1),
        ('north pole', [0,0,0,0,sp.Rational(17,100),-sp.Rational(23,100)],4,1),
        ('cusp interior', [1/mp.sqrt(5)]*4+[.17,-.23],4,1),
    ]
    for label,point in [('phase interpolation',radial_point(mp.mpf('1.07'))),
                        ('south collar',radial_point(mp.mpf('.4'))),
                        ('north collar',radial_point(mp.mpf('2.5')))]:
        values,missing,sign = chart_input(point)
        fixtures.append((label,values,missing,sign))
    graphs = []
    for label, values, missing, sign in fixtures:
        # mpmath receives witness numbers; q0,...,q5 remain free symbols.
        values = [mp.mpf(str(v.evalf(65))) if isinstance(v,sp.Basic) else v for v in values]
        symbolic = D.symbolic_definition(values,missing=missing,sign=sign)
        A, H, assignments = symbolic.differentiate()
        assert A.shape == (6,6) and len(H) == 6
        assert all(h == h.T for h in H)
        assert not any(e.has(sp.Float,sp.Derivative) for _,e in assignments)
        active = symbolic.graph.active_assignments(symbolic.outputs)
        assert not any(n.expression.has(sp.Float) for n in active)
        graphs.append((symbolic,A,H,assignments,values))
        check(f"D5 {label}: actual symbolic graph, {len(active)} values / {len(assignments)} derivative assignments", True)
        with mp.workdps(60):
            symbolic_A,symbolic_H = evaluate_symbolic_jets(symbolic,A,H,assignments,values)
            actual = D.coordinate_tensor(values,missing=missing,sign=sign)
        error_A = np.max(abs(symbolic_A-actual.jacobian))/(1+np.max(abs(actual.jacobian)))
        error_H = np.max(abs(symbolic_H-actual.hessian))/(1+np.max(abs(actual.hessian)))
        check(f"D7 {label}: independent symbolic/numeric derivative comparison (d), {error_A:.2g}/{error_H:.2g}",
              error_A<1e-11 and error_H<1e-10)
    return graphs


def numerical_regions():
    """Safety net (d), never used as a premise of the symbolic identities."""
    points = [radial_point(mp.mpf(r)) for r in
              ('.2','.49','.5','.505','.75','1.07','1.25','1.99','2','5')]
    for inf in (False,True):
        for height in map(mp.mpf,('-.4','-.1','.1','.4')):
            phase = mp.exp(2*mp.pi*mp.j*height)
            points.append(round_point(-mp.conj(phase) if inf else mp.mpc(0),
                                      mp.mpc(0) if inf else phase))
        def cap(X):
            Z = [mp.sqrt(1-X*X)/mp.sqrt(3),
                 mp.sqrt(1-X*X)*mp.sqrt(mp.mpf(2)/3)*mp.exp(mp.j/3)]
            return [Z[0].real,Z[0].imag,Z[1].real,Z[1].imag,X,mp.mpf('.17'),mp.mpf('-.23')]
        for distance in map(mp.mpf,('0','.00000005','.000003','.0000051','.00001')):
            points.append(cap((1-distance)*(1 if inf else -1)))
        for radius in map(mp.mpf,('1.1e-8','1.4e-8')):
            h = (mp.pi if inf else 2*mp.pi)*radius
            points.append(cap((1 if inf else -1)*(1-h)/mp.sqrt(1+h*h)))
    worst_N = worst_square = inverse_error = 0.
    kinds = set()
    formulas = D.NumericalFormulas()
    for point in points:
        u,missing,sign = chart_input(point)
        tensor = D.coordinate_tensor(u,missing=missing,sign=sign)
        kinds.add(tensor.chart)
        worst_N = max(worst_N,np.max(abs(tensor.nijenhuis))/
                      (1+np.max(abs(tensor.matrix))*np.max(abs(tensor.derivative))))
        worst_square = max(worst_square,np.linalg.norm(tensor.matrix@tensor.matrix+np.eye(6))/
                           (1+np.linalg.norm(tensor.matrix)**2))
        assert np.max(abs(tensor.hessian-tensor.hessian.transpose(0,2,1))) < 1e-12*(1+np.max(abs(tensor.hessian)))
        kind,data = formulas.coordinates((point[0]+mp.j*point[1],point[2]+mp.j*point[3]),point[4],point[5:])
        ZZ,XX,phi = formulas.forward(kind,data)
        recovered = [ZZ[0].real,ZZ[0].imag,ZZ[1].real,ZZ[1].imag,XX]
        error = max([abs(x-y) for x,y in zip(point[:5],recovered)]
                    +[abs(mp.exp(2*mp.pi*mp.j*(x-y))-1) for x,y in zip(point[5:],phi)])
        inverse_error = max(inverse_error,float(error))
    check(f"D8 {len(points)} regional tensors and inverse pairs (d); N {worst_N:.3g}, square {worst_square:.3g}, inverse {inverse_error:.3g}",
          len(points)==32 and kinds=={'mid','c0','ci'} and worst_N<1e-9
          and worst_square<1e-11 and inverse_error<1e-11)
    # Five different coordinate computations must yield the same global frame.
    overlap_error = periodic_error = 0.
    for point in ([1/mp.sqrt(5)]*5+[mp.mpf('.17'),mp.mpf('-.23')],
                  radial_point(mp.mpf('1.07')),radial_point(mp.mpf('.4'))):
        reference = D.complex_structure(point)
        for missing in range(5):
            u,_,sign = chart_input(point,missing)
            tensor = D.coordinate_tensor(u,missing=missing,sign=sign)
            T,inverse = D.frame_matrices(point[:5],missing)
            overlap_error = max(overlap_error,np.linalg.norm(T@tensor.matrix@inverse-reference)/(1+np.linalg.norm(reference)))
            J,derivative,frame_N = D.global_frame_data(point[:5],tensor)
            transformed_N = np.einsum('ck,kij,ia,jb->cab',T,tensor.nijenhuis,inverse,inverse)
            scale = 1+np.max(abs(J))*np.max(abs(derivative))+np.max(abs(J))**2
            assert np.max(abs(frame_N-transformed_N))/scale<1e-11
        shifted = point[:5]+[point[5]+2,point[6]-3]
        periodic_error = max(periodic_error,np.linalg.norm(D.complex_structure(shifted)-reference)/(1+np.linalg.norm(reference)))
    for sign in (-1,1):
        J = D.complex_structure([0,0,0,0,sign,.17,-.23])
        assert np.max(abs(J@J+np.eye(6))) < 1e-12
    check(f"D9 all five round charts at three points, torus representatives and public pole calls (d); {overlap_error:.3g}/{periodic_error:.3g}",
          overlap_error<1e-10 and periodic_error<1e-10)
    # The public symbolic interface returns the same global frame convention.
    symbolic = D.complex_structure([0,0,0,0,-1,.17,-.23],symbolic=True)
    assert symbolic.matrix.shape==(6,6) and not symbolic.matrix.has(sp.Float,sp.Derivative)
    # A literal exact evaluation at the 1e-8 gap exceeded a three-second
    # exploratory guard. Numerical tails must return promptly at both signs.
    numerical = D.NumericalFormulas()
    for inverse in (False,True):
        value = mp.mpf('.2')-mp.mpf('1e-8')
        assert numerical._height_R(value,inverse)==value
    # Check the cited VALUE estimate against an independently evaluated tail.
    with mp.workdps(180):
        value = mp.mpf('.199')
        exact = [D.Formulas()._height_R(value,inverse) for inverse in (False,True)]
        bound = 32*mp.mpf('.2')*mp.exp(-mp.mpf('.2')/(mp.mpf('.2')-value))
        with mp.workdps(60):
            approximate = [numerical._height_R(value,inverse) for inverse in (False,True)]
        assert all(abs(x-y)<bound for x,y in zip(exact,approximate))
    check("D11 public symbolic matrix and explicit numerical flat-tail regression (d)",True)


def main():
    start = time.perf_counter()
    mp.mp.dps = 60
    exact_frame()
    exact_tensor()
    exact_phase()
    exact_jet_rules()
    symbolic_primitive_rules()
    formula_correspondence()
    symbolic_cases()
    numerical_regions()
    print(f"ALL PASS ({time.perf_counter()-start:.2f} s)")
    print("Scope: exact identities and the reported source branches; global domain and smoothness contracts require the cited written proofs.")


if __name__ == '__main__':
    main()

