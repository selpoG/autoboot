BeginPackage["ToCpp`"]

cppeval::usage = "cppeval[x]"
createCpp::usage = "createCpp[secs,scalarsecs,vals,numval,eqs,extops,extdelta,filename]"

Begin["`Private`"]

cppeval[n_Integer] := TemplateApply["real(\"``\")", n]
cppeval[a_Times] := StringRiffle[cppeval /@ List @@ a, {"(", " * ", ")"}]
cppeval[a_Plus] := StringRiffle[cppeval /@ List @@ a, {"(", " + ", ")"}]
cppeval[Power[a_, b_]] := TemplateApply["mp::pow(``, ``)", {cppeval[a], cppeval[b]}]
cppeval[Surd[x_, n_]] /; x >= 0 := TemplateApply["mp::pow(``, ``)", {cppeval[x], cppeval[1/n]}]
cppeval[Surd[x_, n_]] := cppeval[-Surd[-x, n]]
cppeval[Sqrt[a_]] := TemplateApply["mp::sqrt(``)", cppeval[a]]
cppeval[Rational[a_, b_]] := TemplateApply["(`` / ``)", {cppeval[a], cppeval[b]}]
cppeval[Sin[a_]] := TemplateApply["mp::sin(``)", cppeval[a]]
cppeval[Cos[a_]] := TemplateApply["mp::cos(``)", cppeval[a]]
cppeval[Tan[a_]] := TemplateApply["mp::tan(``)", cppeval[a]]
cppeval[Sec[a_]] := TemplateApply["mp::sec(``)", cppeval[a]]
cppeval[Csc[a_]] := TemplateApply["mp::csc(``)", cppeval[a]]
cppeval[Cot[a_]] := TemplateApply["mp::cot(``)", cppeval[a]]
cppeval[Exp[a_]] := TemplateApply["mp::exp(``)", cppeval[a]]
cppeval[Log[a_]] := TemplateApply["mp::log(``)", cppeval[a]]
cppeval[Log[b_, a_]] := cppeval[Log[a] / Log[b]]
cppeval[Pi] = "mp::const_pi()";
cppeval[E] = cppeval[Exp[1]];
cppeval[x_?ExactNumberQ] := TemplateApply["real(\"``\")", ToString @ FortranForm @ N[x, 400]]
cppeval[x_?NumberQ] := TemplateApply["real(\"``\")", ToString @ FortranForm[x]]

createCpp[secs_, scalarsecs_, vals_, numval_, eqs_, numext_, extops_, extdelta_, filename_] :=
	TemplateApply[template, <|"secs"->secs, "scalarsecs"->scalarsecs, "vals"->vals, "numval"->numval, "eqs"->eqs, "numext"->numext, "extops"->extops, "extdelta"->extdelta, "filename"->filename|>]

(* filename: "deltas[\"e\"].str(8) + \"-\" + ..." *)
(* numext: 3 *)
(* extdelta: "deltas[\"s\"] = real(args[1]);\n..." *)
(* extops: "ops.emplace(\"s\", Op(deltas[\"s\"], 0, c));\n..." *)
(* scalarsecs: "secs.emplace_back(\"(scalar, 0)\", 2);\n..." *)
(* secs: "{\nSector s(\"even\", 2, ContinuousType);\nfor (const auto& spin: spins) if (spin % 2 == 0) s.add_op(spin);\nsecs.push_back(s);\n}\n..." *)
(* vals: "val[0] = real(-1);\n..." *)
(* eqs: "{\nEq eq(boot, Odd);\neq.add(\"scalar;0\", 0, 0, ops[\"s\"], ext(\"e\", \"s\", \"e\", \"s\"));\neq.add(\"odd+\", ext(\"e\", \"s\", \"e\", \"s\"));\neq.add(\"odd-\", ext(\"e\", \"s\", \"e\", \"s\"));\nboot.add_equation(eq);\n}..." *)

template = "#include <array>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include \"bootstrap_equation.hpp\"
#include \"complex_function.hpp\"
#include \"context.hpp\"
#include \"matrix.hpp\"
#include \"polynomial_program.hpp\"
#include \"primary_op.hpp\"
#include \"real.hpp\"

using algebra::Vector;
using mp::real, mp::rational, mp::parse;
using qboot::Context, qboot::PolynomialProgram, qboot::BootstrapEquation, qboot::Sector;
using std::array, std::map, std::set, std::move, std::vector, std::string, std::unique_ptr;
namespace fs = qboot::fs;

template <class T>
using dict = map<string, T, std::less<>>;
using Op = qboot::PrimaryOperator;
using Eq = qboot::Equation;
constexpr auto ContinuousType = qboot::SectorType::Continuous;
constexpr auto Odd = algebra::FunctionSymmetry::Odd;
constexpr auto Even = algebra::FunctionSymmetry::Even;

static string name(const dict<rational>& deltas);
static PolynomialProgram create(const Context& c, const dict<rational>& deltas, uint32_t numax, set<uint32_t> spins);

string name(const dict<rational>& deltas)
{
	return string(\"sdp-\") + `filename`;
}

PolynomialProgram create(const Context& c, const dict<rational>& deltas, uint32_t numax, set<uint32_t> spins)
{
	dict<Op> ops;
	`extops`
	auto ext = [&ops](auto o1, auto o2, auto o3, auto o4) {
		return array{ops.at(o1), ops.at(o2), ops.at(o3), ops.at(o4)};
	};
	Op one{c};
	vector<Sector> secs;
	// you can add discrete sectors
	// for example, to add sector \"hoge\" whose size of the matrix is 4,
	//   secs.emplace_back(\"hoge\", 4);
	// if you know OPE coefficients of the sector, for example, {1.2, 0.7, -0.1, 0.3},
	//   secs.emplace_back(\"hoge\", 4, Vector{real(\"1.2\"), real(\"0.7\"), real(\"-0.1\"), real(\"0.3\")});
	secs.emplace_back(\"unit\", 1, Vector{real(1)});
	`scalarsecs`
	// you can customize spectrums
	// example:
	//   customize
	//     Sector s(\"hoge\", 2, ContinuousType);
	//     for (const auto& spin: spins)
	//         if (spin % 2 == 0) s.add_op(spin);
	//     secs.push_back(s);
	//   to
	//     Sector s(\"hoge\", 2, ContinuousType);
	//     s.add_op(0, real(\"3\"));    // first scalar: 3 <= delta (irrelevance)
	//     s.add_op(2);               // use unitarity bound for spin-2
	//     s.add_op(4, real(\"8.1\"));  // first spin-4: 8.1 <= delta
	//     // sector \"hoge\" consists of even-spin operators, and we customized spin-0, 2, 4
	//     // use unitarity bound for other operators
	//     for (const auto& spin: spins)
	//         if (spin % 2 == 0 && spin >= 6) s.add_op(spin);
	//     secs.push_back(s);
	`secs`
	// do not edit from here
	Vector<real> val(`numval`);
	`vals`
	BootstrapEquation boot(c, secs, numax);
	`eqs`
	boot.finish();
	// to maximize (resp. minimize) OPE in \"hoge\" sector,
	// call boot.ope_maximize(\"hoge\", \"unit\") (resp. ope_minimize)
	return boot.find_contradiction(\"unit\");
}

int main(int argc, char* argv[])
{
	// internal precision (in binary digits)
	mp::global_prec = 1000;
	mp::global_rnd = MPFR_RNDN;
	// n_Max: the order of taylor expansion of gBlock (we recommend n_Max >= 0.4 * global_prec)
	// lambda: controls the number of derivatives (z = x + sqrt(y), (der x) ^ m (der y) ^ n for m + 2 n <= lambda)
	// dim: the dimension of the physical space
	// numax: controls the number of poles picked in the gBlock (numax + min(spin, numax) / 2)
	constexpr uint32_t n_Max = 400, lambda = 14, dim = 3, numax = 6;
	// spins: spins for the continuous sectors
	set<uint32_t> spins;
	for (uint32_t s = 0; s < 27; ++s) spins.insert(s);
	spins.merge(set{49u, 50u});
	assert(argc == `numext`);
	unique_ptr<char*[]> args(argv);
	dict<rational> deltas;
	// external scalars
	`extdelta`
	args.release();
	Context c(n_Max, lambda, dim);
	auto prob = create(c, deltas, numax, spins);
	auto dir = fs::current_path() / name(deltas);
	move(prob).create_input().write(dir);
	return 0;
}
"

Protect[cppeval, createCpp]

End[ ]

EndPackage[ ]
