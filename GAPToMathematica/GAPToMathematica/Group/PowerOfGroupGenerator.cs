using System.Linq;
using System.Collections.Generic;
using Sprache;

namespace GAPToMathematica
{
	interface IPowerOfGroupGenerator { }
	class PowerOfGroupGenerator : IPowerOfGroupGenerator
	{
		public readonly GroupGenerator Element;
		public readonly int Power;
		public PowerOfGroupGenerator(GroupGenerator x, int n)
		{
			Element = x;
			Power = n;
		}
		public override string ToString() => Power == 1 ? Element.ToString() : $"MatrixPower[{Element}, {Power}]";
		public static readonly Parser<IPowerOfGroupGenerator> Parser;
		static PowerOfGroupGenerator()
		{
			var power = from c in Parse.Char('^')
						from p in Func.PositiveInteger
						select p;
			power = power.Or(x => new Result<int>(1, x));
			Parser = from x in GroupGenerator.Parser
					 from p in power
					 select new PowerOfGroupGenerator(x, p);
		}
	}
	class GroupIdentity : IPowerOfGroupGenerator
	{
		public static readonly GroupIdentity Identity = new GroupIdentity();
		GroupIdentity() { }
		public static readonly Parser<IPowerOfGroupGenerator> Parser = from mark in Parse.String("<identity>of...")
																	   select Identity;
	}
	class GroupElement
	{
		public readonly List<IPowerOfGroupGenerator> Products;
		public GroupElement(IEnumerable<IPowerOfGroupGenerator> p) => Products = new List<IPowerOfGroupGenerator>(p);
		public string ToString(GroupGenerator some) => Products.Count == 0 ? $"MatrixPower[{some}, 0]" : string.Join(".", Products);
		public static readonly Parser<GroupElement> Parser;
		static GroupElement()
		{
			Parser = from x in GroupIdentity.Parser
					 select new GroupElement(Enumerable.Empty<IPowerOfGroupGenerator>());
			Parser = Parser.Or(from x in PowerOfGroupGenerator.Parser.DelimitedBy(Parse.Char('*'))
							   select new GroupElement(x));
		}
	}
}
