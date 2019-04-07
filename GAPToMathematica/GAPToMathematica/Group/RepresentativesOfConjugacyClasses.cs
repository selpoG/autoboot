using System.Linq;
using System.Collections.Generic;
using Sprache;

namespace GAPToMathematica
{
	class RepresentativesOfConjugacyClasses
	{
		public readonly GroupGenerators AllGenerators;
		public readonly List<GroupElement> Representatives;
		public RepresentativesOfConjugacyClasses(GroupGenerators gen, IEnumerable<GroupElement> reps)
		{
			AllGenerators = gen;
			Representatives = new List<GroupElement>(reps);
		}
		public static readonly Parser<RepresentativesOfConjugacyClasses> Parser = from all in GroupGenerators.Parser
																				  from mark in Parse.String("(cggen)=")
																				  from x in GroupElement.Parser.ToSequence()
																				  select new RepresentativesOfConjugacyClasses(all, x);
		public override string ToString() => "{" + string.Join(", ", Representatives.Select(x => x.ToString(AllGenerators.Generators[0]))) + "}";
	}
}
