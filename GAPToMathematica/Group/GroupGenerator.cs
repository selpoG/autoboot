using System.Collections.Generic;
using Sprache;

namespace GAPToMathematica
{
	class GroupGenerator
	{
		static readonly Dictionary<string, GroupGenerator> objs = new Dictionary<string, GroupGenerator>();
		public static GroupGenerator Get(string name)
		{
			if (!objs.ContainsKey(name)) objs[name] = new GroupGenerator(name);
			return objs[name];
		}
		public string Generator;
		GroupGenerator(string name) => Generator = name;
		public static readonly Parser<GroupGenerator> Parser = from a in Func.LowerAlphabet
															   from b in Func.PositiveInteger.Optional()
															   select Get(b.IsDefined ? a + b.Get().ToString() : a + "");
		public override string ToString() => $"\"{Generator}\"";
	}
	class GroupGenerators
	{
		public int Count => Generators.Count;
		public readonly List<GroupGenerator> Generators;
		public GroupGenerators(IEnumerable<GroupGenerator> gen) => Generators = new List<GroupGenerator>(gen);
		public static readonly Parser<GroupGenerators> Parser = from mark in Parse.String("(gen)=")
																from gens in GroupGenerator.Parser.ToSequence()
																select new GroupGenerators(gens);
		public override string ToString() => "{" + string.Join(", ", Generators) + "}";
	}
}
