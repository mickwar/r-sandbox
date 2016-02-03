

nextlevel=c(
    20,     40,     70,     110,    160,
    220,    300,    400,    520,    650,
    790,    940,    1100,   1270,   1450,
    1640,   1840,   2050,   2270,   2500,
    2740,   2990,   3250,   2530,   3800,
    4090,   4390,   4700,   5020,   5350,
    5690,   6040,   6400,   6770,   7150,
    7540,   7940,   8350,   8770,   9200,
    9640,   10090,  10550,  11020,  11500,
    11990,  12490,  13000,  13520,  14050,
    14590,  15140,  15700,  16270,  16850,
    17440,  18040,  18650,  19270,  19900,
    20540,  21190,  21850,  22520,  23200,
    23890,  24590,  25300,  26020,  26750,
    27490,  28240,  29000,  29770,  30550,
    31340,  32140,  32950,  33770,  34600,
    35440,  36290,  37150,  38020,  38900,
    39790,  40690,  41600,  42520,  43450,
    44390,  45340,  46300,  47270,  48250,
    49240,  50240,  51250,  Inf)


update=function(char,addxp,table,status){
    if (status==0)
        return(char)
    if (status==1)
        char[2]=char[2]+floor(addxp*0.75)
    if (status==2)
        char[2]=char[2]+addxp
    multiple=0
    while (multiple==0){
        multiple=1
        temp=char[1]
        if (char[2] >= table[char[1]]){
            char[1]=char[1]+1
            if (status==2)
                char[2]=char[2]-table[char[1]-1]
            else
                char[2]=0
            }
        if (temp != char[1])
            multiple=0
        }
    return(char)
    }

battle=list(
                            #Add: Crono Marle
    list(6,
        c(2,1,0,0,0,0,0)),  #Imp, In: Crono; Not: Marle

                            #Add: Lucca
    list(32,
        c(0,1,2,0,0,0,0),   #Nagette, In: Crono, Lucca;
        c(2,1,0,0,0,0,0)),  #Not: Marle

                            #Add: Frog
    list(53,
        c(2,1,0,0,0,0,0),        #Hench and Diablo, In: Crono, Lucca, Frog
        c(0,1,2,0,0,0,0),        #Not: Marle
        c(0,1,0,2,0,0,0)),
    list(50,    c(2,1,0,0,0,0,0),        #Yakra, In: Crono, Lucca, Frog
            c(0,1,2,0,0,0,0),        #Not: Marle
            c(0,1,0,2,0,0,0)),
    list(40,    c(2,1,0,1,0,0,0),        #Dragon Tank, In: Crono, Lucca
            c(0,1,2,1,0,0,0)),    #Not: Marle, Frog
    list(300,    c(2,0,0,1,0,0,0),        #Gaurdian, In: Crono, Marle, Lucca
            c(0,2,0,1,0,0,0),        #Not: Frog
            c(0,0,2,1,0,0,0)),
    list(54,    c(2,0,0,1,0,0,0),        #Bugger, In: Crono, Marle, Lucca
            c(0,2,0,1,0,0,0),        #Not: Frog
            c(0,0,2,1,0,0,0)),
                            #Add: Robo
    list(33,    c(2,0,1,1,0,0,0),        #Acid
            c(0,2,1,1,0,0,0),        #In: Crono, Robo
            c(0,0,1,1,2,0,0),        #Optional: Marle, Lucca
            c(2,1,0,1,0,0,0),        #Not: Frog
            c(0,1,2,1,0,0,0),
            c(0,1,0,1,2,0,0)),
    list(165,    c(2,0,1,1,0,0,0),        #5 Acids
            c(0,2,1,1,0,0,0),        #In: Crono, Robo
            c(0,0,1,1,2,0,0),        #Optional: Marle, Lucca
            c(2,1,0,1,0,0,0),        #Not: Frog
            c(0,1,2,1,0,0,0),
            c(0,1,0,1,2,0,0)),
    list(480,    c(2,0,1,1,0,0,0),        #R-Series
            c(0,2,1,1,0,0,0),        #In: Crono
            c(2,1,0,1,0,0,0),        #Optional: Marle, Lucca
            c(0,1,2,1,0,0,0)),    #Not: Frog, Dead: Robo?
    list(22,    c(2,0,0,1,1,0,0),        #Hench
            c(0,2,0,1,1,0,0),        #In: Crono
            c(0,0,2,1,1,0,0),        #Optional: Marle, Lucca, Robo
            c(2,0,1,1,0,0,0),        #Not: Frog
            c(0,2,1,1,0,0,0),
            c(0,0,1,1,2,0,0),
            c(2,1,0,1,0,0,0),
            c(0,1,2,1,0,0,0),
            c(0,1,0,1,2,0,0)),
    list(250,    c(2,0,0,1,1,0,0),        #Heckran
            c(0,2,0,1,1,0,0),        #In: Crono
            c(0,0,2,1,1,0,0),        #Optional: Marle, Lucca, Robo
            c(2,0,1,1,0,0,0),        #Not: Frog
            c(0,2,1,1,0,0,0),
            c(0,0,1,1,2,0,0),
            c(2,1,0,1,0,0,0),
            c(0,1,2,1,0,0,0),
            c(0,1,0,1,2,0,0)),
    list(106,    c(2,0,0,1,1,0,0),        #Deceased x2
            c(0,2,0,1,1,0,0),        #In: Crono
            c(0,0,2,1,1,0,0),        #Optional: Marle, Lucca, Robo
            c(2,0,1,1,0,0,0),        #Not: Frog
            c(0,2,1,1,0,0,0),
            c(0,0,1,1,2,0,0),
            c(2,1,0,1,0,0,0),
            c(0,1,2,1,0,0,0),
            c(0,1,0,1,2,0,0)),
    list(144,    c(2,0,0,1,1,0,0),        #Deceased x3
            c(0,2,0,1,1,0,0),        #In: Crono
            c(0,0,2,1,1,0,0),        #Optional: Marle, Lucca, Robo
            c(2,0,1,1,0,0,0),        #Not: Frog
            c(0,2,1,1,0,0,0),
            c(0,0,1,1,2,0,0),
            c(2,1,0,1,0,0,0),
            c(0,1,2,1,0,0,0),
            c(0,1,0,1,2,0,0)),
    list(350,    c(2,0,0,1,1,0,0),        #Zombor
            c(0,2,0,1,1,0,0),        #In: Crono
            c(0,0,2,1,1,0,0),        #Optional: Marle, Lucca, Robo
            c(2,0,1,1,0,0,0),        #Not: Frog
            c(0,2,1,1,0,0,0),
            c(0,0,1,1,2,0,0),
            c(2,1,0,1,0,0,0),
            c(0,1,2,1,0,0,0),
            c(0,1,0,1,2,0,0)),
    list(32,    c(2,0,0,1,1,0,0),        #Goblin
            c(0,2,0,1,1,0,0),        #In: Crono
            c(0,0,2,1,1,0,0),        #Optional: Marle, Lucca, Robo
            c(2,0,1,1,0,0,0),        #Not: Frog
            c(0,2,1,1,0,0,0),
            c(0,0,1,1,2,0,0),        #Consecutive Battle
            c(2,1,0,1,0,0,0),
            c(0,1,2,1,0,0,0),
            c(0,1,0,1,2,0,0)),
    list(32,    c(2,0,0,1,1,0,0),        #Ogan
            c(0,2,0,1,1,0,0),        #In: Crono
            c(0,0,2,1,1,0,0),        #Optional: Marle, Lucca, Robo
            c(2,0,1,1,0,0,0),        #Not: Frog
            c(0,2,1,1,0,0,0),
            c(0,0,1,1,2,0,0),        #Consecutive Battle
            c(2,1,0,1,0,0,0),
            c(0,1,2,1,0,0,0),
            c(0,1,0,1,2,0,0)),
    list(400,    c(2,0,0,1,1,0,0),        #Masa & Mune
            c(0,2,0,1,1,0,0),        #In: Crono
            c(0,0,2,1,1,0,0),        #Optional: Marle, Lucca, Robo
            c(2,0,1,1,0,0,0),        #Not: Frog
            c(0,2,1,1,0,0,0),
            c(0,0,1,1,2,0,0),        #Consecutive Battle
            c(2,1,0,1,0,0,0),
            c(0,1,2,1,0,0,0),
            c(0,1,0,1,2,0,0)),
    list(500,    c(2,0,0,1,1,0,0),        #Masamune
            c(0,2,0,1,1,0,0),        #In: Crono
            c(0,0,2,1,1,0,0),        #Optional: Marle, Lucca, Robo
            c(2,0,1,1,0,0,0),        #Not: Frog
            c(0,2,1,1,0,0,0),
            c(0,0,1,1,2,0,0),        #Consecutive Battle
            c(2,1,0,1,0,0,0),
            c(0,1,2,1,0,0,0),
            c(0,1,0,1,2,0,0)),
    list(44,    c(2,0,0,1,1,0,0),        #Nereides
            c(0,2,0,1,1,0,0),        #In: Crono
            c(0,0,2,1,1,0,0),        #Optional: Marle, Lucca, Robo
            c(2,0,1,1,0,0,0),        #Not: Frog
            c(0,2,1,1,0,0,0),
            c(0,0,1,1,2,0,0),
            c(2,1,0,1,0,0,0),
            c(0,1,2,1,0,0,0),
            c(0,1,0,1,2,0,0)),
    list(360,    c(2,0,0,1,1,0,0),        #Reptites
            c(0,2,0,1,1,0,0),        #In: Crono
            c(0,0,2,1,1,0,0),        #Optional: Marle, Lucca, Robo
            c(2,0,1,1,0,0,0),        #Not: Frog
            c(0,2,1,1,0,0,0),
            c(0,0,1,1,2,0,0),        #Consecutive Battle
            c(2,1,0,1,0,0,0),
            c(0,1,2,1,0,0,0),
            c(0,1,0,1,2,0,0)),
    list(288,    c(2,0,0,1,1,0,0),        #Reptites
            c(0,2,0,1,1,0,0),        #In: Crono
            c(0,0,2,1,1,0,0),        #Optional: Marle, Lucca, Robo
            c(2,0,1,1,0,0,0),        #Not: Frog
            c(0,2,1,1,0,0,0),
            c(0,0,1,1,2,0,0),        #Consecutive Battle
            c(2,1,0,1,0,0,0),
            c(0,1,2,1,0,0,0),
            c(0,1,0,1,2,0,0)),
                            #Add: Ayla
    list(162,    c(2,0,1,1,1,0,0),        #Evilweevil
            c(0,2,1,1,1,0,0),        #In: Crono, Ayla
            c(0,0,1,1,1,2,0),        #Optional: Marle, Lucca, Robo
            c(2,1,0,1,1,0,0),        #Not: Frog
            c(0,1,2,1,1,0,0),
            c(0,1,0,1,1,2,0),
            c(2,1,1,1,0,0,0),
            c(0,1,1,1,2,0,0),
            c(0,1,1,1,0,2,0)),
    list(147,    c(2,0,1,1,1,0,0),        #Evilweevil
            c(0,2,1,1,1,0,0),        #In: Crono, Ayla
            c(0,0,1,1,1,2,0),        #Optional: Marle, Lucca, Robo
            c(2,1,0,1,1,0,0),        #Not: Frog
            c(0,1,2,1,1,0,0),
            c(0,1,0,1,1,2,0),
            c(2,1,1,1,0,0,0),
            c(0,1,1,1,2,0,0),
            c(0,1,1,1,0,2,0)),
    list(500,    c(2,0,1,1,1,0,0),        #Nizbel
            c(0,2,1,1,1,0,0),        #In: Crono, Ayla
            c(0,0,1,1,1,2,0),        #Optional: Marle, Lucca, Robo
            c(2,1,0,1,1,0,0),        #Not: Frog
            c(0,1,2,1,1,0,0),
            c(0,1,0,1,1,2,0),
            c(2,1,1,1,0,0,0),
            c(0,1,1,1,2,0,0),
            c(0,1,1,1,0,2,0)),
    list(488,    c(2,0,1,0,1,1,0),        #Hench + Vamp
            c(0,2,1,0,1,1,0),        #In: Crono, Frog
            c(0,0,1,2,1,1,0),        #Optional: Marle, Lucca, Robo
            c(2,1,0,0,1,1,0),        #Not: Ayla
            c(0,1,2,0,1,1,0),
            c(0,1,0,2,1,1,0),
            c(2,1,1,0,0,1,0),
            c(0,1,1,0,2,1,0),
            c(0,1,1,2,0,1,0)),
    list(60,    c(2,0,1,0,1,1,0),        #Decedant
            c(0,2,1,0,1,1,0),        #In: Crono, Frog
            c(0,0,1,2,1,1,0),        #Optional: Marle, Lucca, Robo
            c(2,1,0,0,1,1,0),        #Not: Ayla
            c(0,1,2,0,1,1,0),
            c(0,1,0,2,1,1,0),
            c(2,1,1,0,0,1,0),
            c(0,1,1,0,2,1,0),
            c(0,1,1,2,0,1,0)),
    list(500,    c(2,0,1,0,1,1,0),        #Slash
            c(0,2,1,0,1,1,0),        #In: Crono, Frog
            c(0,0,1,2,1,1,0),        #Optional: Marle, Lucca, Robo
            c(2,1,0,0,1,1,0),        #Not: Ayla
            c(0,1,2,0,1,1,0),
            c(0,1,0,2,1,1,0),
            c(2,1,1,0,0,1,0),
            c(0,1,1,0,2,1,0),
            c(0,1,1,2,0,1,0)),
    list(500,    c(2,0,1,0,1,1,0),        #Flea
            c(0,2,1,0,1,1,0),        #In: Crono, Frog
            c(0,0,1,2,1,1,0),        #Optional: Marle, Lucca, Robo
            c(2,1,0,0,1,1,0),        #Not: Ayla
            c(0,1,2,0,1,1,0),
            c(0,1,0,2,1,1,0),
            c(2,1,1,0,0,1,0),
            c(0,1,1,0,2,1,0),
            c(0,1,1,2,0,1,0)),
    list(72,    c(2,0,1,0,1,1,0),        #Decedant
            c(0,2,1,0,1,1,0),        #In: Crono, Frog
            c(0,0,1,2,1,1,0),        #Optional: Marle, Lucca, Robo
            c(2,1,0,0,1,1,0),        #Not: Ayla
            c(0,1,2,0,1,1,0),
            c(0,1,0,2,1,1,0),
            c(2,1,1,0,0,1,0),
            c(0,1,1,0,2,1,0),
            c(0,1,1,2,0,1,0)),
    list(462,    c(2,0,1,0,1,1,0),        #Outlaw + Groupie
            c(0,2,1,0,1,1,0),        #In: Crono, Frog
            c(0,0,1,2,1,1,0),        #Optional: Marle, Lucca, Robo
            c(2,1,0,0,1,1,0),        #Not: Ayla
            c(0,1,2,0,1,1,0),
            c(0,1,0,2,1,1,0),
            c(2,1,1,0,0,1,0),
            c(0,1,1,0,2,1,0),
            c(0,1,1,2,0,1,0)),
    list(512,    c(2,0,1,0,1,1,0),        #Juggler
            c(0,2,1,0,1,1,0),        #In: Crono, Frog
            c(0,0,1,2,1,1,0),        #Optional: Marle, Lucca, Robo
            c(2,1,0,0,1,1,0),        #Not: Ayla
            c(0,1,2,0,1,1,0),
            c(0,1,0,2,1,1,0),
            c(2,1,1,0,0,1,0),
            c(0,1,1,0,2,1,0),
            c(0,1,1,2,0,1,0)),
    list(1500,    c(2,0,1,0,1,1,0),        #Magus
            c(0,2,1,0,1,1,0),        #In: Crono, Frog
            c(0,0,1,2,1,1,0),        #Optional: Marle, Lucca, Robo
            c(2,1,0,0,1,1,0),        #Not: Ayla
            c(0,1,2,0,1,1,0),
            c(0,1,0,2,1,1,0),
            c(2,1,1,0,0,1,0),
            c(0,1,1,0,2,1,0),
            c(0,1,1,2,0,1,0)),
    list(880,    c(2,0,1,1,1,0,0),        #Nizbel II
            c(0,2,1,1,1,0,0),        #In: Crono, Ayla
            c(0,0,1,1,1,2,0),        #Optional: Marle, Lucca, Frog, Robo
            c(2,1,0,1,1,0,0),    
            c(0,1,2,1,1,0,0),
            c(0,1,0,1,1,2,0),
            c(2,1,1,0,1,0,0),
            c(0,1,1,2,1,0,0),
            c(0,1,1,0,1,2,0),
            c(2,1,1,1,0,0,0),
            c(0,1,1,1,2,0,0),
            c(0,1,1,1,0,2,0)),
    list(1800,    c(2,0,1,1,1,0,0),        #Azala + Black Tyrano
            c(0,2,1,1,1,0,0),        #In: Crono, Ayla
            c(0,0,1,1,1,2,0),        #Optional: Marle, Lucca, Frog, Robo
            c(2,1,0,1,1,0,0),    
            c(0,1,2,1,1,0,0),
            c(0,1,0,1,1,2,0),
            c(2,1,1,0,1,0,0),
            c(0,1,1,2,1,0,0),
            c(0,1,1,0,1,2,0),
            c(2,1,1,1,0,0,0),
            c(0,1,1,1,2,0,0),
            c(0,1,1,1,0,2,0))
    )

com=function(x,battle){
    x=x-1
    combos=rep(0,length(battle))
    multiples=rep(0,length(battle))
    multiples[1]=1
    for (i in 2:length(battle)){
        multiples[i]=multiples[i-1]*(length(battle[[i-1]])-1)
        }
    multiples[1]=0
    for (j in length(battle):1){
        if (j>1){
            if (multiples[j] != multiples[j-1]){
                if (x>=multiples[j]){
                    combos[j]=floor(x/multiples[j])
                    x=x-multiples[j]*combos[j]
                    }
                }
            }
        }
    return(combos+1)
    }

firstcombo=integer(35)+1

nextcombo=function(combo,battle){
    combo[1]=combo[1]+1
    for (j in 2:length(battle)){
        if (combo[j-1] > length(battle[[j-1]])-1){
            combo[j-1]=1
            combo[j]=combo[j]+1
            }
        }
    if (combo[35] > length(battle[35]))
        print("Maxed out")
    return(combo)
    }

(combo=firstcombo)
while (1){
    (combo=nextcombo(combo,battle))
    print(combo)
    }
### Order is: Crono, Marle, Lucca, Frog, Robo, Ayla, Magus
characters=list(c(1,0),c(1,0),c(2,0),c(5,0),c(10,0),c(18,0),c(37,0))

addxp=function(characters,addxp,table,status){
    for (i in 1:7){
        characters[[i]]=update(characters[[i]],addxp,table,status[i])
        }
    return(characters)
    }

### Simulation

# first 10000 combos
norders = 10000
endlevels=matrix(0,norders,8)

endlevels=matrix(0,norders,8)

# up to how many battles used
upto = length(battle)

### keepers
keep = list(double(length(battle)), double(8))

# batnum = death combo number
for (batnum in 1:norders){
    characters=list(c(1,0),c(1,0),c(2,0),c(5,0),c(10,0),c(18,0),c(37,0))
    if (batnum==1){
        #endlevels=matrix(0,1,8)
        order=com(1,battle)
        }
    if (batnum>1){
        order=com(batnum,battle)
        }
    for (i in 1:upto){
        characters=addxp(characters,battle[[i]][[1]],
            nextlevel,battle[[i]][[order[i]+1]])
        }
    for (k in 1:7){
        endlevels[batnum,k]=characters[[k]][1]
        }
    endlevels[batnum,8]=sum(endlevels[batnum,1:7])
    }

which.min(endlevels[,8])
endlevels[1,]

com(4962300000000000000000000001, battle)

nextcombo(com(4962300000000000000000000003, battle), battle)

### Total number of battle combinations
ugh = prod(unlist(lapply(battle, function(x) length(x)-1)))

ugh / 60 / 60 / 24 / 365 / 1000000000000
