source("./card_functions.R")

mickey = make_deck("")
mickey = add_card(mickey, "+1")
mickey = add_card(mickey, "+1")
mickey = add_card(mickey, "+1")
mickey = add_card(mickey, "+1")
mickey = add_card(mickey, "+1")
mickey = add_card(mickey, "pierce 3;redraw")
mickey = add_card(mickey, "pierce 3;redraw")
mickey = add_card(mickey, "disarm;redraw")
mickey = add_card(mickey, "muddle;redraw")

jesse = make_deck("")
jesse = rm_card(jesse, "-1")
jesse = rm_card(jesse, "-1")
jesse = rm_card(jesse, "-1")
jesse = rm_card(jesse, "-1")
jesse = rm_card(jesse, "-1")
jesse = rm_card(jesse, "-2")
jesse = add_card(jesse, "+0")
jesse = add_card(jesse, "+1")
jesse = add_card(jesse, "+1;redraw")
jesse = add_card(jesse, "+1;redraw")

enemy = make_deck("")
ooze = make_deck("ooze")
drake = make_deck("rending drake")
living_spirit = make_deck("living spirit")


# Modifier
mickey = draw_card(mickey)
#mickey = shuffle_deck(mickey)

jesse = draw_card(jesse)
#jesse = shuffle_deck(jesse)

enemy = draw_card(enemy)
#enemy = shuffle_deck(enemy)

# Turns
drake = draw_card(drake, TRUE)
living_spirit = draw_card(living_spirit, TRUE)
ooze = draw_card(ooze, TRUE)

