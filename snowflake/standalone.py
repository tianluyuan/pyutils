def simple_hese_names(run_id, event_id):
    names = {
        (118178,66452255):'Camilla the Chicken',	#3
        (118283,9445773):'Rowlf',	#5
        (118381,19162840):'Beaker',	#7
        (118435,58198553):'Mr. Snuffleupagus',	#9
        (118545,63733662):'Bert',	#11
        (118549,11722208):'Sweetums',	#13
        (118602,23096391):'Dr. Bunsen Honeydew',	#15
        (118607,40435683):'Lew Zealand',	#17
        (119196,37713300):'Zoot',	#19
        (119214,8606380):'Link Hogthrob',	#21
        (119316,36556705):'Ernie',	#23
        (119352,56498321):'Guy Smiley',	#25
        (119404,80750561):'Miss Piggy',	#27
        (119470,48424887):'Gonzo the Great',	#29
        (119474,33152537):'Rizzo',	#31
        (119595,30769232):'Count von Count',	#33
        (119674,8449256):'Kermit',	#35
        (119842,82622124):'Sam Eagle',	#37
        (120045,22615214):'Dr. Teeth',	#39
        (115994,2538090):'Fozzy',	#41
        (115994,29874216):'Scooter',	#43
        (116528,52433389):'Animal',	#45
        (116698,10198436):'Swedish Chef',	#47
        (116876,63208734):'Dr. Strangepork',	#49
        (116878,34919596):'Statler and Waldorf',	#51
        (117322,7422546):'Pepe the King Prawn',	#53
        (117371,31623515):'Crazy Harry',	#55
        (117782,49441871):'Floyd',	#57
        (118145,5142726):'Oscar',	#59
        (121240,72944671):'Big Bird'	
    }
    try:
        return f'{run_id} {event_id} ({names[run_id, event_id]})'
    except KeyError:
        return f'{run_id} {event_id}'
