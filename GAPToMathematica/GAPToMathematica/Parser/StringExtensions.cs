using System;
using System.Collections.Generic;
using System.Linq;

namespace Sprache
{
	static class StringExtensions
	{
		public static IEnumerable<char> ToEnumerable(this string x)
		{
#if STRING_IS_ENUMERABLE
			return x;
#else
			if (x == null) throw new ArgumentNullException("this");
			for (var i = 0; i < x.Length; ++i) yield return x[i];
#endif
		}

		public static string Join<T>(string separator, IEnumerable<T> values)
		{
#if STRING_JOIN_ENUMERABLE
            return string.Join(separator, values);
#else
			return string.Join(separator, values.Select(v => v.ToString()).ToArray());
#endif
		}
	}
}
