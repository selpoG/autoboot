using System.Collections.Generic;
using System.Linq;
using Sprache;

namespace GAPToMathematica
{
	class SizeOfConjugacyClasses
	{
		public int Count;
		public int[] Size;
		public SizeOfConjugacyClasses(IEnumerable<int> cgs)
		{
			Size = cgs.ToArray();
			Count = Size.Length;
		}
		public static readonly Parser<SizeOfConjugacyClasses> Parser = from mark in Parse.String("(cgsize)=")
																	   from x in Func.PositiveInteger.ToSequence()
																	   select new SizeOfConjugacyClasses(x);
		public override string ToString() => "{" + string.Join(", ", Size) + "}";
	}
}
