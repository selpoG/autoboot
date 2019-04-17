using System;
using System.Linq;
using System.Collections.Generic;
using Sprache;

namespace GAPToMathematica
{
	class GroupGenerators
	{
		public int Count => Generators.Count;
		public readonly List<string> Generators;
		public GroupGenerators(IEnumerable<FreeGroupElement> gen)
		{
			Generators = new List<string>();
			foreach (var g in gen)
			{
				if (g.Products.Count != 1) throw new ArgumentException($"{g} is not a primitive.");
				foreach (var p in g.Products)
				{
					if (p.p != 1) Console.Error.WriteLine($"Warning: x{p.x}^({p.p}) was given as a generator. Please check x{p.x} is involutive.");
					Generators.Add($"x{p.x}");
				}
			}
		}
		public static readonly Parser<GroupGenerators> Parser = from mark in Parse.String("(gen)=")
																from gens in FreeGroupElement.Parser.ToSequence()
																select new GroupGenerators(gens);
		public override string ToString() => "{" + string.Join(", ", Generators.Select(s => $"\"{s}\"")) + "}";
	}
}
