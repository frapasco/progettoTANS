# tans

Come girare il programma (cosa compilare, ordine dei comandi, parametri da passare a funzione/mettere in macro); parametri di default
cosa e come fa; scopo; Fisica (riassunto); parametri e variabili principali; spiegazione delle formule usate
illustrare classi ( e macro)
tree
risultati (con graficci) e conclusioni
commenti







NOTA: per la ricostruzione si considerano solo eventi con almeno due tracce intere vere dentro il rivelatore: così, in presenza di noise anche a angoli phi tali da poter generare una tracklet (deltaphi fra due punti di noise (uno per rivelatore)<phimax), è possibile trovare un accumulo di intersezioni delle tracklet nei dintorni di ztrue con attorno picchetti spuri dal noise (fatto salvo eventi particolari).
Per come è fatto l'salgoritmo, nel caso di due accumuli analoghi (pari contenuto del bin dell'istogramma) viene selezionato il rpimo, che può accidentamente corrispondere a ztrue anche se in realtà il vertice non è ben ricostruito (bisogna scegli9ere a caso uno dei due bin e sperare sia giusto), quindi l'efficenza viene in parte sovrastimata. Ma non ricostruendo eventi a una traccia dentro il rivelatore, oltre a eliminare eventi con n bin tutti con contenuto 1 (noise a phi tali da generrare tracklet) si eliminano anche evneti in cui la posizione di hit e noise è tale che viene costruita una sola tracklet (corrispondente alla traccia vera) e quindi il vertice potrebbe essere ricavato: si ha quinddi una riduzione dell'efficenza (si perdono degli eventi buoni).
La almeno parziale compensazione fra i due fattori (aumento e riduzione efficenza) porterà presumibilmente a avvicinarsi di più all'efficenza vera rispetto al considerare solo uno dei due fattori 




Nota: nella simulazione si è scelto di considerare (salvandone gli hit) solo le tracce che geometricamente intersecano entrambi i layer: a livello di riciostruzione questo ha poca influenza in quanto le tracce che intersecano solo un layer comportano punti su uno strato per cui non esiste nessun punto sull’altro strato taloe che una retta fra essi punti al vertice e si comportano quindi per la ricostruzione come i punti di noise che sono aggiunti (da noi con distribuzione uniforme ma il costruttore intpint(x,y,z) permetterebbe di metterli in posizioni a piacere) nella macro ricostruzione. Si alleggerisce così a livello computazionale il calcolo di phimax (1 ciclo for invece di due innestati). Nella classe intpoint è comunque presente il dm label che si può usare per salvare il numero di traccia a cui appartiene ogni singolo intpoint e quindi le macro utente sono modificabili per considerare nella simulazione anche le tracce a un solo hit. Vengono inoltre ignorati (e.g. nella realtà con un trigger) i verticei tali che nessun prodotto passa entrami i riverlatori(e quindi eventi inutili per il vertexing)  e quindi l’efficienza trovata non conterrà la componente di efficienza legata ai fqattori geometrici e cinematici di evento e apparato che comportano il non avere nessun prodotto cfhe passa entrami i rivelatori (ma essendo le efficienzae m9oltiplicative si può calcolare questa componente di efficienza e fare il prodotto se interressati)
