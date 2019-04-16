using System.Collections.Generic;
using System.Linq;
using Sprache;

namespace GAPToMathematica
{
	// ±a/b*E(n)^m => (±a, b, n, m)
	// (a,b) = 1, b > 0
	// n > 1 => 0 < m < n
	// n = 1 => m = 0
	// Parser assumes that a > 0 (recognize '+', '-' as delimiters)
	// "5", "1/4", "E(6)", "E(3)^2", "2*E(6)", "3*E(4)^3", "1/3*E(5)", "5/7*E(8)^5"
	class AlgebraicTerm
	{
		public int Numerator { get; private set; }
		public readonly int Denomitor, DegreeOfRoot, Power;
		AlgebraicTerm(int n, int d, AlgebraicTerm rp) : this(n, d)
		{
			if (rp == null)
			{
				DegreeOfRoot = 1;
				Power = 0;
			}
			else
			{
				DegreeOfRoot = rp.DegreeOfRoot;
				Power = rp.Power;
			}
		}
		public AlgebraicTerm(int n, int d)
		{
			Numerator = n;
			Denomitor = d;
			DegreeOfRoot = 1;
			Power = 0;
		}
		public AlgebraicTerm(int n, int d, int r, int p) : this(n, d)
		{
			DegreeOfRoot = r;
			Power = p;
		}
		public AlgebraicTerm Negate()
		{
			Numerator *= -1;
			return this;
		}
		public override string ToString()
		{
			var num = Denomitor == 1 ? Numerator.ToString() : $"{Numerator}/{Denomitor}";
			if (DegreeOfRoot == 1) return num;
			var exp = Power == 1 ? $"e[{DegreeOfRoot}]" : $"e[{DegreeOfRoot}, {Power}]";
			switch (num)
			{
				case "1": return exp;
				case "-1": return "-" + exp;
				default: return $"{num}*{exp}";
			}
		}
		public string ToStringWithSign()
		{
			var x = ToString();
			return x[0] == '-' ? x : '+' + x;
		}
		public static readonly Parser<AlgebraicTerm> Parser;
		static AlgebraicTerm()
		{
			Parser = from e in Parse.String("E(")
					 from r in Func.PositiveInteger
					 from f in Parse.Char(')')
					 from p in (from m in Parse.Char('^')
								from p in Func.PositiveInteger
								select p).Optional()
					 select new AlgebraicTerm(1, 1, r, p.IsDefined ? p.Get() : 1);
			Parser = Parser.Or(from n in Func.PositiveInteger
							   from d in (from m in Parse.Char('/')
										  from d in Func.PositiveInteger
										  select d).Optional()
							   from rp in (from t in Parse.Char('*')
										   from rp in Parser
										   select rp).Optional()
							   select new AlgebraicTerm(n, d.IsDefined ? d.Get() : 1, rp.GetOrDefault()));
		}
	}
	class AlgebraicExpression
	{
		public readonly List<AlgebraicTerm> Terms;
		public AlgebraicExpression(IEnumerable<AlgebraicTerm> ts) => Terms = new List<AlgebraicTerm>(ts);
		public void Add(AlgebraicTerm x) => Terms.Add(x);
		public override string ToString()
		{
			switch (Terms.Count)
			{
				case 0: return "0";
				case 1: return Terms[0].ToString();
				default: return Terms[0] + string.Join("", Terms.Skip(1).Select(x => x.ToStringWithSign()));
			}
		}
		public static readonly Parser<AlgebraicExpression> Parser;
		static AlgebraicExpression()
		{
			Parser = from st in Parse.Char('-').Optional()
					 from x in AlgebraicTerm.Parser
					 let nx = st.IsDefined ? x.Negate() : x
					 from y in (from s in Func.Sign
								from z in AlgebraicTerm.Parser
								select s == '+' ? z : z.Negate()).Many()
					 select new AlgebraicExpression(x.Concat(y));
			Parser = Parser.Or(from zero in Parse.Char('0')
							   select new AlgebraicExpression(Enumerable.Empty<AlgebraicTerm>()));
		}
	}
}
