using System.Linq;
using System.Collections.Generic;
using Sprache;
using System.Diagnostics;

namespace GAPToMathematica
{
	class FreeGroupElement
	{
		readonly List<(int x, int p)> products;
		public IReadOnlyList<(int x, int p)> Products => products;
		public FreeGroupElement(IList<int> l)
		{
			Debug.Assert(l.Count % 2 == 0);
			products = new List<(int x, int p)>(l.Count / 2);
			for (var i = 0; i < l.Count; i += 2) products.Add((l[i], l[i + 1]));
		}
		public override string ToString() => products.Count == 0 ? $"MatrixPower[x1, {0}]" : string.Join(".", products.Select(p => $"MatrixPower[x{p.x}, {p.p}]"));
		public static readonly Parser<FreeGroupElement> Parser = from x in Func.ToSequence(Func.Integer)
																 select new FreeGroupElement(x);
	}
}
