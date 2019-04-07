using System.Text;
using System.Linq;
using Sprache;

namespace GAPToMathematica
{
	// non-pcgroup
	// [60, 5], [120, 5], [120, 34], [120, 35], [168, 42], [180, 19], [240, 89], [240, 90], [240, 91], [240, 92], [240, 93], [240, 94], [240, 189], [240, 190],
	// [300, 22], [336, 114], [336, 208], [336, 209], [360, 51], [360, 118], [360, 119], [360, 120], [360, 121], [360, 122], [420, 13],
	// [480, 217], [480, 218], [480, 219], [480, 220], [480, 221], [480, 222], [480, 943], [480, 944], [480, 945], [480, 946], [480, 947], [480, 948], [480, 949], [480, 950],
	// [480, 951], [480, 952], [480, 953], [480, 954], [480, 955], [480, 956], [480, 957], [480, 958], [480, 959], [480, 960], [480, 1186], [480, 1187],
	// for i in [1..m] do if not IsPcGroup(SmallGroup(n,i)) then Print("[", n, ", ", i, "]"); fi; od;
	partial class Group
	{
		public readonly GroupID Id;
		public readonly SizeOfConjugacyClasses Sizes;
		public readonly CharacterTable Characters;
		public readonly IrreducibleRepresentations Irreps;
		public readonly RepresentativesOfConjugacyClasses Representatives;
		public readonly GroupGenerators Generators;
		public readonly int NumberOfConjugacyClasses, NumberOfGenerators;
		public Group(GroupID id, SizeOfConjugacyClasses scg, CharacterTable ct, IrreducibleRepresentations irrep)
		{
			Id = id;
			Sizes = scg;
			Characters = ct;
			Irreps = irrep;
			Representatives = Irreps.Representatives;
			Generators = Representatives.AllGenerators;
			NumberOfConjugacyClasses = Sizes.Count;
			NumberOfGenerators = Generators.Count;
		}
		public bool Check()
		{
			if (Id.Order < 2 || Id.Id < 1) return false;
			if (NumberOfConjugacyClasses < 1 || Sizes.Size.Any(x => x <= 0)) return false;
			if (Characters.Characters.Count != NumberOfConjugacyClasses) return false;
			if (Characters.Characters.Any(x => x.CharacterElements.Count != NumberOfConjugacyClasses)) return false;
			if (NumberOfGenerators < 1) return false;
			if (Representatives.Representatives.Count != NumberOfConjugacyClasses) return false;
			if (Irreps.Irreps.Count != NumberOfConjugacyClasses) return false;
			if (Irreps.Irreps.Any(x => x.Count != NumberOfGenerators)) return false;
			foreach (var ir in Irreps.Irreps)
			{
				var d = ir[0].GetLength(0);
				if (d < 1) return false;
				foreach (var m in ir) if (d != m.GetLength(0)) return false;
			}
			return true;
		}
		public override string ToString()
		{
			var sb = new StringBuilder();
			sb.Append($"(id) = {Id}\n");
			sb.Append($"(cgsize) = {Sizes}\n");
			sb.Append($"(chartab) = {Characters}\n");
			sb.Append($"(gen) = {Generators}\n");
			sb.Append($"(cggen) = {Representatives}\n");
			sb.Append($"(irrep) = {Irreps}");
			return sb.ToString();
		}
		public static readonly Parser<Group> Parser = from id in GroupID.Parser
													  from cgs in SizeOfConjugacyClasses.Parser
													  from ct in CharacterTable.Parser
													  from ir in IrreducibleRepresentations.Parser
													  select new Group(id, cgs, ct, ir);
	}
}
