add_card = function(deck, card){
    deck = rbind(deck, c(card, "0"))
    return (deck)
    }
rm_card = function(deck, card){
    ind = which(deck[,1] == card)
    if (length(ind) == 0)
        stop("Card not found in deck: ", card)
    deck = deck[-min(ind),]  # min() to account for multiple copies
    return (deck)
    }
draw_card = function(deck, autoshuffle = FALSE){
    draw_possible = which(deck[,2] == "0")
    if (length(draw_possible) == 1){
        draw_ind = draw_possible
    } else {
        draw_ind = sample(draw_possible, 1)
        }
    card_text = strsplit(deck[draw_ind,1], ";")[[1]]

    # If initiative
    has_init = grep("init", card_text)
    init = ""
    if (length(has_init) > 0){
        #init = strsplit(card_text[has_init], " ")[[1]][2]
        init = card_text[has_init]
        card_text = card_text[-has_init]
        }

    # text of card
    text = card_text[1]
    card_text = card_text[-1]

    # If shuffle or redraw
    special = ""
    if (length(card_text) > 0)
        special = card_text[1]

    deck[draw_ind, 2] = "1" # put card in discard pile
    print(c(init, text, special))
    if (special == "redraw")
        return (draw_card(deck))    # recursive if redraw

    if (autoshuffle && special == "shuffle")
        deck[,2] = "0"

    return (deck)
    }

# Not part of draw_card function since shuffle_deck occurs at the end
# the round, not immediately when a shuffle card is drawn
shuffle_deck = function(deck){
    deck[,2] = "0"
    return (deck)
    }

make_deck = function(type  = ""){
    type = tolower(type)
    type = gsub(" ", "_", type)
    # Beginning player/monster modification deck
    if (type == ""){
        x = NULL
        x = add_card(x, "2x;shuffle")
        x = add_card(x, "+2")
        for (i in 1:5)
            x = add_card(x, "+1")
        for (i in 1:6)
            x = add_card(x, "+0")
        for (i in 1:5)
            x = add_card(x, "-1")
        x = add_card(x, "-2")
        x = add_card(x, "miss;shuffle")
    } else {
        list.files("./enemy/")
        paste0("./enemy/", type, ".txt")
        card_lines = as.character(read.table(paste0("./enemy/", type, ".txt"), sep = "$", quote = "$")[,1])
        x = NULL
        for (i in 1:length(card_lines))
            x = add_card(x, card_lines[i])
        }
    if (missing(x))
        stop("type not found")
    return (x)
    }
