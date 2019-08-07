//Arjun Singh Gambhir
//Code ported from latscat, with many changes included as comments.
#ifndef _MAPSTUFF
#define _MAPSTUFF

//map with default value:
template <template<class,class,class...> class C, typename K, typename V, typename... Args>
V GetWithDef(const C<K,V,Args...>& m, K const& key, const V& defval){
	typename C<K,V,Args...>::const_iterator it = m.find( key );
	if (it == m.end()) return defval;
	return it->second;
}

#endif
