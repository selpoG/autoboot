using System.Linq;
using System.Collections.Generic;
using Sprache;

namespace GAPToMathematica
{
	class RepresentativesOfConjugacyClasses
	{
		public readonly GroupGenerators AllGenerators;
		public readonly List<FreeGroupElement> Representatives;
		public RepresentativesOfConjugacyClasses(GroupGenerators gen, IEnumerable<FreeGroupElement> reps)
		{
			AllGenerators = gen;
			Representatives = new List<FreeGroupElement>(reps);
		}
		public static readonly Parser<RepresentativesOfConjugacyClasses> Parser = from all in GroupGenerators.Parser
																				  from mark in Parse.String("(cggen)=")
																				  from x in FreeGroupElement.Parser.ToSequence()
																				  select new RepresentativesOfConjugacyClasses(all, x);
		public override string ToString() => "{" + string.Join(", ", Representatives) + "}";
	}
}
