## code to prepare `species200` dataset goes here

# this is an example dataset, including a vector of 200 species from the Ross
# Sea (Antarctica, Southern Ocean), that yield NCBI records in the numbers
# between 10 and 30 (approximately), with fasta sequences never exceeding 2000
# base pairs (bp) in length, and never below 100 bp (at 12 Jun 2024).

species200 <- c("Aetideopsis minor", "Newnesia antarctica", "Ophiacantha densispina",
                "Gaetanus tenuispinus", "Fungiacyathus marenzelleri", "Cinachyra antarctica",
                "Pallenopsis vanhoeffeni", "Caryophyllia diomedeae", "Tremaster mirabilis",
                "Euphausia triacantha", "Trophon shackletoni", "Labidiaster annulatus",
                "Paracucumis turricata", "Prosipho spiralis", "Trypanedenta gigantea",
                "Abatus nimrodi", "Puncturella spirigera", "Gorekia crassicirris",
                "Vanadis antarctica", "Liothyrella neozelanica", "Mesothuria bifurcata",
                "Fannyella spinosa", "Proserpinaster neozelanicus", "Epimeria robusta",
                "Charcotia obesa", "Trichobranchus glacialis", "Macroscapha turbida",
                "Metridia gerlachei", "Eunice pennata", "Haliclona scotti", "Crucella scotiae",
                "Protelpidia murrayi", "Notisis elongata", "Nymphon charcoti",
                "Fenestrulina malusii", "Ophiura rouchi", "Ammothea longispina",
                "Orchomenella franklini", "Scaphocalanus brevicornis", "Themisto gaudichaudii",
                "Cucumaria georgiana", "Odontaster meridionalis", "Laevilacunaria antarctica",
                "Amphicteis gunneri", "Graneledone antarctica", "Acodontaster conspicuus",
                "Sclerasterias mollis", "Polymastia invaginata", "Bathydraco scotiae",
                "Starkus mixtus", "Amphiura belgicae", "Macroptychaster accrescens",
                "Edentoliparis terraenovae", "Leucon antarcticus", "Antarctotetilla sagitta",
                "Heterocucumis denticulata", "Molpadia musculus", "Harmothoe acuminata",
                "Psalidaster mordax", "Altenaeum charcoti", "Scolecithricella minor",
                "Ophiocten megaloplax", "Amphiura joubini", "Paradiplospinus gracilis",
                "Amythas membranifera", "Austrolaenilla antarctica", "Calanoides acutus",
                "Epimeria inermis", "Placiphorella atlantica", "Phyllodoce longipes",
                "Aricidea belgicae", "Cinachyra barbata", "Primnoisis formosa",
                "Antarctodomus thielei", "Onogorgia nodosa", "Convexella magelhaenica",
                "Novocrania huttoni", "Atolla wyvillei", "Rossella racovitzae",
                "Cryothenia amphitreta", "Rhincalanus gigas", "Megaleledone setebos",
                "Heterocucumis steineni", "Probuccinum tenerum", "Antarctotetilla leptoderma",
                "Parastenella spinosa", "Rossella antarctica", "Iothia emarginuloides",
                "Thouarella variabilis", "Tentorium papillatum", "Cirratulus cirratus",
                "Tritoniella belli", "Trachythyone bouvetensis", "Philobrya wandelensis",
                "Staurocucumis liouvillei", "Pasiphaea acutifrons", "Antarctinoe ferox",
                "Enypniastes eximia", "Spiophanes kroyeri", "Clausocalanus laticeps",
                "Ophiocamax gigas", "Paraliparis antarcticus", "Thysanoessa macrura",
                "Gymnoscopelus hintonoides", "Propeamussium investigatoris",
                "Mastigoteuthis psychrophila", "Pogonophryne macropogon", "Aforia magnifica",
                "Gymnoscopelus piabilis", "Bathydraco macrolepis", "Torellia insignis",
                "Bathylagus antarcticus", "Acanthonotozomoides oatesi", "Notomastus latericeus",
                "Spinocalanus magnus", "Diplasterias brandti", "Amphitrite kerguelensis",
                "Perknaster densus", "Cucamba psolidiformis", "Pachycara brachycephalum",
                "Adamussium colbecki", "Poecillastra compressa", "Pseudorchomene rossi",
                "Stauroteuthis gilchristi", "Oithona frigida", "Rhachotropis abyssalis",
                "Epimeria schiaparelli", "Antarctothoa bougainvillei", "Doto antarctica",
                "Thouarella pendulina", "Chlanidota signeyana", "Primnoella antarctica",
                "Pogonophryne cerebropogon", "Sterechinus antarcticus", "Hymenodora gracilis",
                "Ammothea spinosa", "Ophioplinthus anceps", "Torellia mirabilis",
                "Silicula rouchi", "Psilaster charcoti", "Cuenotaster involutus",
                "Rossella nuda", "Notobdella nototheniae", "Corella eumyota",
                "Eulagisca gigantea", "Cynomacrurus piriei", "Echinopsolus charcoti",
                "Abyssocucumis abyssorum", "Laonice cirrata", "Ophiuroglypha irrorata",
                "Echinopsolus mollis", "Cyamiomactra laminifera", "Neobuccinum eatoni",
                "Maldane sarsi", "Aglaophamus trissophyllus", "Peniagone vignoni",
                "Colossendeis tortipalpis", "Pasiphaea scotiae", "Astrochlamys bruneus",
                "Antarctophiline alata", "Gymnoscopelus opisthopterus", "Limatula hodgsoni",
                "Pareledone panchroma", "Tentorium semisuberites", "Alternatipathes alternata",
                "Nacella kerguelenensis", "Antarctotetilla grandis", "Primnoisis gracilis",
                "Sphaerotylus capitatus", "Hippomedon kergueleni", "Sphaerotylus antarcticus",
                "Spinocalanus abyssalis", "Lebbeus kiae", "Eunoe opalina", "Arntzia gracilis",
                "Gymnoscopelus bolini", "Zoroaster spinulosus", "Molpadiodemas villosus",
                "Protomyctophum bolini", "Triconia conifera", "Abatus cavernosus",
                "Dasystenella acanthina", "Synoicum adareanum", "Antimargarita dulcis",
                "Muraenolepis microps", "Epimeria rimicarinata", "Pogonophryne marmorata",
                "Psolicrux coatsi", "Latrunculia biformis", "Acodontaster hodgsoni",
                "Paraeuchaeta antarctica", "Halipteris willemoesi", "Eulagisca uschakovi",
                "Lysasterias perrieri", "Rhopiella hirsuta", "Ophiocten dubium",
                "Paramphinome australis", "Hyperiella dilatata", "Abatus shackletoni",
                "Hyperiella macronyx", "Epimeria macronyx", "Waegelea antarctica",
                "Colossendeis scotti", "Bathyplotes bongraini", "Spinocalanus antarcticus",
                "Ophiosteira echinulata", "Ophiosparte gigas")

usethis::use_data(species200, overwrite = TRUE)
