using Sprache;

namespace GAPToMathematica
{
	class GroupID
	{
		public int Order, Id;
		public GroupID(int g, int i) { Order = g; Id = i; }
		public static readonly Parser<GroupID> Parser = from mark in Parse.String("(id)=")
														from x in Func.PositiveInteger.ToSequence()
														select new GroupID(x[0], x[1]);
		public override string ToString() => "{" + $"{Order}, {Id}" + "}";
	}
}
