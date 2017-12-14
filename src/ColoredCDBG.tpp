void Color::join(const Color& data, CompactedDBG<Color>& cdbg){

    ColoredCDBG& colored_cdbg = static_cast<ColoredCDBG&>(cdbg);
}

void Color::split(const size_t pos, const size_t len, Color& new_data, CompactedDBG<Color>& cdbg) const {

    new_data.color_id = color_id;
}
