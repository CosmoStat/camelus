void add_bias{cosmo_hm *cmhm,gal_map *gMap, gal_map *gMap_bias, error **err)
{
	int i,j;
	double b=0.00856 ;

	for (i=0; i<gMap->length; i++) {
		gList = gMap->map[i];
		for (j=0, gNode=gList->first; j<gList->size; j++, gNode=gNode->next) {
			g=gNode->g;
			

}
