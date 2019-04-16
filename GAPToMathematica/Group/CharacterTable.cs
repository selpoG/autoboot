using System.Collections.Generic;
using System.Linq;
using Sprache;

namespace GAPToMathematica
{
	class Character
	{
		public readonly List<AlgebraicExpression> CharacterElements;
		public Character(IEnumerable<AlgebraicExpression> z) => CharacterElements = new List<AlgebraicExpression>(z);
		public override string ToString() => "{" + string.Join(", ", CharacterElements) + "}";
		public static readonly Parser<Character> Parser = from x in AlgebraicExpression.Parser.ToSequence()
														  select new Character(x);
	}
	class CharacterTable
	{
		public readonly List<Character> Characters;
		public CharacterTable(IEnumerable<Character> cs) => Characters = new List<Character>(cs);
		public override string ToString() => "{" + string.Join(", ", Characters) + "}";
		public static readonly Parser<CharacterTable> Parser;
		static CharacterTable() => Parser = from mark in Parse.String("(chartab)=")
											from cs in (from c in Character.Parser
														select c).AtLeastOnce()
											select new CharacterTable(cs);
	}
}
