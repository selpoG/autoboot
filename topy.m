BeginPackage["ToPython`"]

pyeval::usage = "pyeval[x]"
createPython::usage = "createPython[secs,scalarnum,vals,rmats,smats,umats,filename]"

Begin["`Private`"]

pyeval[n_Integer] := TemplateApply["context(\"``\")", n]
pyeval[a_Times] := StringRiffle[pyeval /@ List @@ a, {"(", " * ", ")"}]
pyeval[a_Plus] := StringRiffle[pyeval /@ List @@ a, {"(", " + ", ")"}]
pyeval[Power[a_, b_]] := TemplateApply["pow(``, ``)", {pyeval[a], pyeval[b]}]
pyeval[Surd[x_, n_]] /; x >= 0 := TemplateApply["pow(``, ``)", {pyeval[x], pyeval[1/n]}]
pyeval[Surd[x_, n_]] := pyeval[-Surd[-x, n]]
pyeval[Sqrt[a_]] := TemplateApply["sqrt(``)", pyeval[a]]
pyeval[Rational[a_, b_]] := TemplateApply["(`` / ``)", {pyeval[a], pyeval[b]}]
pyeval[Pi] = "context(pi)";
pyeval[E] = "context(e)";
pyeval[EulerGamma] = "context(euler_gamma)";
pyeval[Catalan] = "context(catalan)";
pyeval[Khinchin] = "context(khinchin)";
pyeval[Glaisher] = "context(glaisher)";
pyeval[Sin[a_]] := TemplateApply["sin(``)", pyeval[a]]
pyeval[Cos[a_]] := TemplateApply["cos(``)", pyeval[a]]
pyeval[Tan[a_]] := TemplateApply["tan(``)", pyeval[a]]
pyeval[Sec[a_]] := TemplateApply["sec(``)", pyeval[a]]
pyeval[Csc[a_]] := TemplateApply["csc(``)", pyeval[a]]
pyeval[Cot[a_]] := TemplateApply["cot(``)", pyeval[a]]
pyeval[Sinh[a_]] := TemplateApply["sinh(``)", pyeval[a]]
pyeval[Cosh[a_]] := TemplateApply["cosh(``)", pyeval[a]]
pyeval[Tanh[a_]] := TemplateApply["tanh(``)", pyeval[a]]
pyeval[Sech[a_]] := TemplateApply["sech(``)", pyeval[a]]
pyeval[Csch[a_]] := TemplateApply["csch(``)", pyeval[a]]
pyeval[Coth[a_]] := TemplateApply["coth(``)", pyeval[a]]
pyeval[ArcSin[a_]] := TemplateApply["asin(``)", pyeval[a]]
pyeval[ArcCos[a_]] := TemplateApply["acos(``)", pyeval[a]]
pyeval[ArcTan[a_]] := TemplateApply["atan(``)", pyeval[a]]
pyeval[ArcSec[a_]] := TemplateApply["asec(``)", pyeval[a]]
pyeval[ArcCsc[a_]] := TemplateApply["acsc(``)", pyeval[a]]
pyeval[ArcCot[a_]] := TemplateApply["acot(``)", pyeval[a]]
pyeval[ArcSinh[a_]] := TemplateApply["asinh(``)", pyeval[a]]
pyeval[ArcCosh[a_]] := TemplateApply["acosh(``)", pyeval[a]]
pyeval[ArcTanh[a_]] := TemplateApply["atanh(``)", pyeval[a]]
pyeval[ArcSech[a_]] := TemplateApply["asech(``)", pyeval[a]]
pyeval[ArcCsch[a_]] := TemplateApply["acsch(``)", pyeval[a]]
pyeval[ArcCoth[a_]] := TemplateApply["acoth(``)", pyeval[a]]
pyeval[Exp[a_]] := TemplateApply["exp(``)", pyeval[a]]
pyeval[Log[a_]] := TemplateApply["log(``)", pyeval[a]]
pyeval[Log[b_, a_]] := TemplateApply["log(``, ``)", {pyeval[a], pyeval[b]}]
pyeval[I] = "I";
pyeval[x_?ExactNumberQ] := TemplateApply["context(\"``\")", Echo[FortranForm @ N[Echo[x], 300]]]
pyeval[x_?NumberQ] := TemplateApply["context(\"``\")", FortranForm[x]]

createPython[secs_, scalarnum_, vals_, rmats_, smats_, umats_, filename_] :=
	TemplateApply[template, <|"secs"->secs, "scalarnum"->scalarnum, "vals"->vals, "rmats"->rmats, "smats"->smats, "umats"->umats, "filename"->filename|>]

template = "# -*- coding: utf-8 -*-
from __future__ import print_function
import sage.cboot as cb
from sage.misc.cachefunc import cached_function
import os.path
import sys
from sage.rings.rational import Rational
from sage.all import pi, e, euler_gamma, catalan, khinchin, glaisher, sin, cos, tan, sec, csc, cot, sinh, cosh, tanh, sech, csch, coth, asin, acos, atan, asec, acsc, acot, asinh, acosh, atanh, asech, acsch, acoth, sqrt, log, exp, I

context = cb.context_for_scalar(epsilon=0.5, Lambda=11)
spins = list(range(22))
secs = {`secs`}
scalarnum = `scalarnum`
nu_max = 8
mygap = {}
def gaps(deltas):
	return mygap

save_dir = \".\"

def get_path(path):
	return os.path.join(save_dir, path)

@cached_function
def prepare_g(spin, Delta_12, Delta_34):
	return context.approx_cb(nu_max, spin, Delta_1_2=Delta_12, Delta_3_4=Delta_34)

def get_shift(gaps, sector, spin):
	if (sector, spin) in gaps: return context(gaps[(sector, spin)])
	elif spin == 0: return context.epsilon
	else: return 2 * context.epsilon + spin

val = `vals`
one = context(1)

def make_F_range(delta, sector, num, spin):
	shift = get_shift(gaps(delta), sector, spin)
	sign = one if spin % 2 == 0 else -one
	F = lambda d1, d2, d3, d4: sign * context.dot(context.F_minus_matrix((d2 + d3) / 2), prepare_g(spin, d1 - d2, d3 - d4).shift(shift))
	H = lambda d1, d2, d3, d4: sign * context.dot(context.F_plus_matrix((d2 + d3) / 2), prepare_g(spin, d1 - d2, d3 - d4).shift(shift))
	get = lambda fh, d1, d2, d3, d4: fh(delta[d1], delta[d2], delta[d3], delta[d4])
`rmats`
	raise RuntimeError(\"unknown sector name\")

def make_F_scalar(delta, num):
	F = lambda d1, d2, d3, d4, d: context.dot(context.F_minus_matrix((d2 + d3) / 2), context.gBlock(0, d, d1 - d2, d3 - d4))
	H = lambda d1, d2, d3, d4, d: context.dot(context.F_plus_matrix((d2 + d3) / 2), context.gBlock(0, d, d1 - d2, d3 - d4))
	get = lambda fh, d1, d2, d3, d4, d: fh(delta[d1], delta[d2], delta[d3], delta[d4], delta[d])
`smats`
	raise RuntimeError(\"unknown sector number\")

def make_F_unit(delta):
	F = lambda d1, d2, d3, d4: context.dot(context.F_minus_matrix((d2 + d3) / 2), context.gBlock(0, 0, d1 - d2, d3 - d4))
	H = lambda d1, d2, d3, d4: context.dot(context.F_plus_matrix((d2 + d3) / 2), context.gBlock(0, 0, d1 - d2, d3 - d4))
	get = lambda fh, d1, d2, d3, d4: fh(delta[d1], delta[d2], delta[d3], delta[d4])
	`umats`

def make_SDP(delta):
	cdel = dict()
	for k in delta: cdel[k] = context(delta[k])
	pvms = []
	for sec in secs:
		for spin in spins:
			if spin % 2 == sec[1]:
				for num in range(secs[sec]): pvms.append(make_F_range(cdel, sec[0], num, spin))
	for num in range(scalarnum): pvms.append(make_F_scalar(cdel, num))
	norm = make_F_unit(cdel)
	obj = 0
	return context.sumrule_to_SDP(norm, obj, pvms)

def name(deltas):
	return \"sdp-`filename`\".format(deltas).replace(\"/\", \"#\")

def has_done(deltas):
	return os.path.exists(get_path(name(deltas) + \".xml\"))

def writefl(mes):
	print(mes, end=\"\")
	sys.stdout.flush()

def write_SDP(deltas):
	prob = get_path(name(deltas) + \".xml\")
	if not has_done(deltas): make_SDP(deltas).write(prob)
"
End[ ]

EndPackage[ ]
