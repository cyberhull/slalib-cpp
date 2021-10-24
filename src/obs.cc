/*
 * C++ Port of the SLALIB library.
 * Written by Vadim Sytnikov.
 * Copyright (C) 2021 CyberHULL, Ltd.
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 */
#include "slalib.h"
#include <cstring>

namespace sla {

/**
 * Retrieves parameters of selected ground-based observing stations.
 *
 * Programs can obtain a list of all currently supported stations by calling the function repeatedly, with
 * `n` = 0,1,2,3... When `false` is returned from the function, the list of stations has been exhausted.
 *
 * Station numbers, identifiers, names and other details are subject to change and should not be hardwired into
 * application programs.
 *
 * IMPORTANT -- BEWARE OF THE LONGITUDE SIGN CONVENTION. The longitude returned by sla::obs() is west-positive in
 * accordance with astronomical usage. However, this sign convention is left-handed and is the opposite of the one
 * used by geographers; elsewhere in SLALIB the preferable east-positive convention is used. In particular, note
 * that for use in sla::aop(), sla::aoppa() and sla::oap() the sign of the longitude must be reversed.
 *
 * @param n Number specifying observing station; if `n` is less than 0, it is ignored, and observatory selection is
 *   performed based on it's `id` (second argument); if `n` is greater than the number of observatories known to the
 *   function, then `false` is returned.
 * @param id Identifier specifying observing station; must be upper-case; if partial `id` is provided, then the first
 *   observatory with full `id` starting with provided partial `id` will be returned; the `id` argument is *only*
 *   considered if supplied `n` is negative.
 * @param obs Return value: parameters of specified observing station.
 * @return `true` if parameters were successfully retrieved, `false` otherwise; in the latter case, `o_id` and
 *   `o_name` fields of the `obs` structure are set to "?" strings, and all numeric fields to 0.0.
 */
bool obs(int n, const char* id, Observatory& obs) {

    auto west = [](int degrees, int arcmins, double arcsecs) constexpr -> double {
        constexpr double ARCSECS_2_RADIANS = 0.484813681109535994e-5;
        return ARCSECS_2_RADIANS * ((double) (60 * (60 * degrees + arcmins)) + arcsecs);
    };
    auto north = [&west](int degrees, int arcmins, double arcsecs) constexpr -> double {
        return west(degrees, arcmins, arcsecs);
    };
    auto east = [&west](int degrees, int arcmins, double arcsecs) constexpr -> double {
        return -west(degrees,arcmins,arcsecs);
    };
    auto south = [&west](int degrees, int arcmins, double arcsecs) constexpr -> double {
        return -west(degrees,arcmins,arcsecs);
    };

    constexpr int N_OBS = 85;
    static const Observatory observatories[N_OBS] = {
        // AAT (Observer's Guide)
        {"AAT",        "Anglo-Australian 3.9m Telescope",         east(149, 3, 57.91),   south(31, 16, 37.34),   1164.0},
        // WHT (Gemini, April 1987)
        {"LPO4.2",     "William Herschel 4.2m Telescope",         west(17, 52, 53.9),    north(28, 45, 38.1),    2332.0},
        // INT (Gemini, April 1987)
        {"LPO2.5",     "Isaac Newton 2.5m Telescope",             west(17, 52, 39.5),    north(28, 45, 43.2),    2336.0},
        // JKT (Gemini, April 1987)
        {"LPO1",       "Jacobus Kapteyn 1m Telescope",            west(17, 52, 41.2),    north(28, 45, 39.9),    2364.0},
        // Lick 120" (S.L.Allen, private communication, 2002)
        {"LICK120",    "Lick 120 inch",                           west(121, 38, 13.689), north(37, 20, 34.931),  1286.0},
        // MMT 6.5m conversion (MMT Observatory website)
        {"MMT",        "MMT 6.5m, Mt Hopkins",                    west(110, 53, 4.4),    north(31, 41, 19.6),    2608.0},
        // Victoria B.C. 1.85m (1984 Almanac)
        {"DAO72",      "DAO Victoria BC 1.85 metre",              west(123, 25, 1.18),   north(48, 31, 11.9),    238.0},
        // Las Campanas (1983 Almanac)
        {"DUPONT",     "Du Pont 2.5m Telescope, Las Campanas",    west(70, 42, 9.0),     south(29, 0, 11.0),     2280.0},
        // Mt Hopkins 1.5m (1983 Almanac)
        {"MTHOP1.5",   "Mt Hopkins 1.5 metre",                    west(110, 52, 39.0),   north(31, 40, 51.4),    2344.0},
        // Mt Stromlo 74" (1983 Almanac)
        {"STROMLO74",  "Mount Stromlo 74 inch",                   east(149, 0, 27.59),   south(35, 19, 14.3),    767.0},
        // ANU 2.3m, SSO (Gary Hovey)
        {"ANU2.3",     "Siding Spring 2.3 metre",                 east(149, 3, 40.3),    south(31, 16, 24.1),    1149.0},
        // Greenbank 140' (1983 Almanac)
        {"GBVA140",    "Greenbank 140 foot",                      west(79, 50, 9.61),    north(38, 26, 15.4),    881.0},
        // Cerro Tololo 4m (1982 Almanac)
        {"TOLOLO4M",   "Cerro Tololo 4 metre",                    west(70, 48, 53.6),    south(30, 9, 57.8),     2235.0},
        // Cerro Tololo 1.5m (1982 Almanac)
        {"TOLOLO1.5M", "Cerro Tololo 1.5 metre",                  west(70, 48, 54.5),    south(30, 9, 56.3),     2225.0},
        // Tidbinbilla 64m (1982 Almanac)
        {"TIDBINBLA",  "Tidbinbilla 64 metre",                    east(148, 58, 48.2),   south(35, 24, 14.3),    670.0},
        // Bloemfontein 1.52m (1981 Almanac)
        {"BLOEMF",     "Bloemfontein 1.52 metre",                 east(26, 24, 18.0),    south(29, 2, 18.0),     1387.0},
        // Bosque Alegre 1.54m (1981 Almanac)
        {"BOSQALEGRE", "Bosque Alegre 1.54 metre",                west(64, 32, 48.0),    south(31, 35, 53.0),    1250.0},
        // USNO 61" astrographic reflector, Flagstaff (1981 Almanac)
        {"FLAGSTF61",  "USNO 61 inch astrograph, Flagstaff",      west(111, 44, 23.6),   north(35, 11, 2.5),     2316.0},
        // Lowell 72" (1981 Almanac)
        {"LOWELL72",   "Perkins 72 inch, Lowell",                 west(111, 32, 9.3),    north(35, 5, 48.6),     2198.0},
        // Harvard 1.55m (1981 Almanac)
        {"HARVARD",    "Harvard College Observatory 1.55m",       west(71, 33, 29.32),   north(42, 30, 19.0),    185.0},
        // Okayama 1.88m (1981 Almanac)
        {"OKAYAMA",    "Okayama 1.88 metre",                      east(133, 35, 47.29),  north(34, 34, 26.1),    372.0},
        // Kitt Peak Mayall 4m (1981 Almanac)
        {"KPNO158",    "Kitt Peak 158 inch",                      west(111, 35, 57.61),  north(31, 57, 50.3),    2120.0},
        // Kitt Peak 90 inch (1981 Almanac)
        {"KPNO90",     "Kitt Peak 90 inch",                       west(111, 35, 58.24),  north(31, 57, 46.9),    2071.0},
        // Kitt Peak 84 inch (1981 Almanac)
        {"KPNO84",     "Kitt Peak 84 inch",                       west(111, 35, 51.56),  north(31, 57, 29.2),    2096.0},
        // Kitt Peak 36 foot (1981 Almanac)
        {"KPNO36FT",   "Kitt Peak 36 foot",                       west(111, 36, 51.12),  north(31, 57, 12.1),    1939.0},
        // Kottamia 74" (1981 Almanac)
        {"KOTTAMIA",   "Kottamia 74 inch",                        east(31, 49, 30.0),    north(29, 55, 54.0),    476.0},
        // La Silla 3.6m (1981 Almanac)
        {"ESO3.6",     "ESO 3.6 metre",                           west(70, 43, 36.0),    south(29, 15, 36.0),    2428.0},
        // Mauna Kea 88 inch (IfA website, Richard Wainscoat)
        {"MAUNAK88",   "Mauna Kea 88 inch",                       west(155, 28, 9.96),   north(19, 49, 22.77),   4213.6},
        // UKIRT (IfA website, Richard Wainscoat)
        {"UKIRT",      "UK Infra Red Telescope",                  west(155, 28, 13.18),  north(19, 49, 20.75),   4198.5},
        // Quebec 1.6m (1981 Almanac)
        {"QUEBEC1.6",  "Quebec 1.6 metre",                        west(71, 9, 9.7),      north(45, 27, 20.6),    1114.0},
        // Mt Ekar 1.82m (1981 Almanac)
        {"MTEKAR",     "Mt Ekar 1.82 metre",                      east(11, 34, 15.0),    north(45, 50, 48.0),    1365.0},
        // Mt Lemmon 60" (1981 Almanac)
        {"MTLEMMON60", "Mt Lemmon 60 inch",                       west(110, 42, 16.9),   north(32, 26, 33.9),    2790.0},
        // Mt Locke 2.7m (1981 Almanac)
        {"MCDONLD2.7", "McDonald 2.7 metre",                      west(104, 1, 17.60),   north(30, 40, 17.7),    2075.0},
        // Mt Locke 2.1m (1981 Almanac)
        {"MCDONLD2.1", "McDonald 2.1 metre",                      west(104, 1, 20.1),    north(30, 40, 17.7),    2075.0},
        // Palomar 200" (1981 Almanac)
        {"PALOMAR200", "Palomar 200 inch",                        west(116, 51, 50.0),   north(33, 21, 22.0),    1706.0},
        // Palomar 60" (1981 Almanac)
        {"PALOMAR60",  "Palomar 60 inch",                         west(116, 51, 31.0),   north(33, 20, 56.0),    1706.0},
        // David Dunlap 74" (1981 Almanac)
        {"DUNLAP74",   "David Dunlap 74 inch",                    west(79, 25, 20.0),    north(43, 51, 46.0),    244.0},
        // Haute Provence 1.93m (1981 Almanac)
        {"HPROV1.93",  "Haute Provence 1.93 metre",               east(5, 42, 46.75),    north(43, 55, 53.3),    665.0},
        // Haute Provence 1.52m (1981 Almanac)
        {"HPROV1.52",  "Haute Provence 1.52 metre",               east(5, 42, 43.82),    north(43, 56, 0.2),     667.0},
        // San Pedro Martir 83" (1981 Almanac)
        {"SANPM83",    "San Pedro Martir 83 inch",                west(115, 27, 47.0),   north(31, 2, 38.0),     2830.0},
        // Sutherland 74" (1981 Almanac)
        {"SAAO74",     "Sutherland 74 inch",                      east(20, 48, 44.3),    south(32, 22, 43.4),    1771.0},
        // Tautenburg 2m (1981 Almanac)
        {"TAUTNBG",    "Tautenburg 2 metre",                      east(11, 42, 45.),     north(50, 58, 51.),     331.0},
        // Catalina 61" (1981 Almanac)
        {"CATALINA61", "Catalina 61 inch",                        west(110, 43, 55.1),   north(32, 25, 0.7),     2510.0},
        // Steward 90" (1981 Almanac)
        {"STEWARD90",  "Steward 90 inch",                         west(111, 35, 58.24),  north(31, 57, 46.9),    2071.0},
        // Russian 6m (1981 Almanac)
        {"USSR6",      "USSR 6 metre",                            east(41, 26, 30.0),    north(43, 39, 12.0),    2100.0},
        // Arecibo 1000' (1981 Almanac)
        {"ARECIBO",    "Arecibo 1000 foot",                       west(66, 45, 11.1),    north(18, 20, 36.6),    496.0},
        // Cambridge 5km (1981 Almanac)
        {"CAMB5KM",    "Cambridge 5km",                           east(0, 2, 37.23),     north(52, 10, 12.2),    17.0},
        // Cambridge 1 mile (1981 Almanac)
        {"CAMB1MILE",  "Cambridge 1 mile",                        east(0, 2, 21.64),     north(52, 9, 47.3),     17.0},
        // Bonn 100m (1981 Almanac)
        {"EFFELSBERG", "Effelsberg 100 metre",                    east(6, 53, 1.5),      north(50, 31, 28.6),    366.0},
        // Greenbank 300' (1981 Almanac) [R.I.P.]
        {"GBVA300",    "Greenbank 300 foot",                      west(79, 50, 56.36),   north(38, 25, 46.3),    894.0},
        // Jodrell Bank Mk 1 (1981 Almanac)
        {"JODRELL1",   "Jodrell Bank 250 foot",                   west(2, 18, 25.),      north(53, 14, 10.5),    78.0},
        // Australia Telescope Parkes Observatory (Peter te Lintel Hekkert)
        {"PARKES",     "Parkes 64 metre",                         east(148, 15, 44.3591),south(32, 59, 59.8657), 391.79},
        // VLA (1981 Almanac)
        {"VLA",        "Very Large Array",                        west(107, 37, 3.82),   north(34, 4, 43.5),     2124.0},
        // Sugar Grove 150' (1981 Almanac)
        {"SUGARGROVE", "Sugar Grove 150 foot",                    west(79, 16, 23.0),    north(38, 31, 14.0),    705.0},
        // Russian 600' (1981 Almanac)
        {"USSR600",    "USSR 600 foot",                           east(41, 35, 25.5),    north(43, 49, 32.0),    973.0},
        // Nobeyama 45 metre mm dish (based on 1981 Almanac entry)
        {"NOBEYAMA",   "Nobeyama 45 metre",                       east(138, 29, 12.0),   north(35, 56, 19.0),    1350.0},
        // James Clerk Maxwell 15 metre mm telescope, Mauna Kea [From GPS measurements on 11Apr2007 for eSMA setup (R.Tilanus)]
        {"JCMT",       "JCMT 15 metre",                           west(155, 28, 37.3),   north(19, 49, 22.22),   4124.75},
        // ESO 3.5 metre NTT, La Silla (K.Wirenstrand)
        {"ESONTT",     "ESO 3.5 metre NTT",                       west(70, 43, 7.0),     south(29, 15, 30.0),    2377.0},
        // St Andrews University Observatory (1982 Almanac)
        {"ST.ANDREWS", "St Andrews",                              west(2, 48, 52.5),     north(56, 20, 12.0),    30.0},
        // Apache Point 3.5 metre (R.Owen)
        {"APO3.5",     "Apache Point 3.5m",                       west(105, 49, 11.56),  north(32, 46, 48.96),   2809.0},
        // W.M.Keck Observatory, Telescope 1 (William Lupton)
        {"KECK1",      "Keck 10m Telescope #1",                   west(155, 28, 28.99),  north(19, 49, 33.41),   4160.0},
        // Tautenberg Schmidt (1983 Almanac)
        {"TAUTSCHM",   "Tautenberg 1.34 metre Schmidt",           east(11, 42, 45.0),    north(50, 58, 51.0),    331.0},
        // Palomar Schmidt (1981 Almanac)
        {"PALOMAR48",  "Palomar 48-inch Schmidt",                 west(116, 51, 32.0),   north(33, 21, 26.0),    1706.0},
        // UK Schmidt, Siding Spring (1983 Almanac)
        {"UKST",       "UK 1.2 metre Schmidt, Siding Spring",     east(149, 4, 12.8),    south(31, 16, 27.8),    1145.0},
        // Kiso Schmidt, Japan (1981 Almanac)
        {"KISO",       "Kiso 1.05 metre Schmidt, Japan",          east(137, 37, 42.2),   north(35, 47, 38.7),    1130.0},
        // ESO Schmidt, La Silla (1981 Almanac)
        {"ESOSCHM",    "ESO 1 metre Schmidt, La Silla",           west(70, 43, 46.5),    south(29, 15, 25.8),    2347.0},
        // Australia Telescope Compact Array (WGS84 coordinates of Station 35, Mark Calabretta)
        {"ATCA",       "Australia Telescope Compact Array",       east(149, 33, 0.5),    south(30, 18, 46.385),  236.9},
        // Australia Telescope Mopra Observatory (Peter te Lintel Hekkert)
        {"MOPRA",      "ATNF Mopra Observatory",                  east(149, 5, 58.732),  south(31, 16, 4.451),   850.0},
        // Subaru telescope, Mauna Kea (IfA website, Richard Wainscoat)
        {"SUBARU",     "Subaru 8m telescope",                     west(155, 28, 33.67),  north(19, 49, 31.81),   4163.0},
        // Canada-France-Hawaii Telescope, Mauna Kea (IfA website, Richard Wainscoat)
        {"CFHT",       "Canada-France-Hawaii 3.6m Telescope",     west(155, 28, 7.95),   north(19, 49, 30.91),   4204.1},
        // W.M.Keck Observatory, Telescope 2 (William Lupton)
        {"KECK2",      "Keck 10m Telescope #2",                   west(155, 28, 27.24),  north(19, 49, 35.62),   4159.6},
        // Gemini North, Mauna Kea (IfA website, Richard Wainscoat)
        {"GEMININ",    "Gemini North 8-m telescope",              west(155, 28, 8.57),   north(19, 49, 25.69),   4213.4},
        // Five College Radio Astronomy Observatory (Tim Jenness)
        {"FCRAO",      "Five College Radio Astronomy Obs",        west(72, 20, 42.0),    north(42, 23, 30.0),    314.0},
        // NASA Infra Red Telescope Facility (IfA website, Richard Wainscoat)
        {"IRTF",       "NASA IR Telescope Facility, Mauna Kea",   west(155, 28, 19.2),   north(19, 49, 34.39),   4168.1},
        // Caltech Submillimeter Observatory (IfA website, Richard Wainscoat; height estimated)
        {"CSO",        "Caltech Sub-mm Observatory, Mauna Kea",   west(155, 28, 31.79),  north(19, 49, 20.78),   4080.0},
        // ESO VLT, UT1 (ESO website, VLT Whitebook Chapter 2)
        {"VLT1",       "ESO VLT, Paranal, Chile: UT1",            west(70, 24, 11.642),  south(24, 37, 33.117),  2635.43},
        // ESO VLT, UT2 (ESO website, VLT Whitebook Chapter 2)
        {"VLT2",       "ESO VLT, Paranal, Chile: UT2",            west(70, 24, 10.855),  south(24, 37, 31.465),  2635.43},
        // ESO VLT, UT3 (ESO website, VLT Whitebook Chapter 2)
        {"VLT3",       "ESO VLT, Paranal, Chile: UT3",            west(70, 24, 9.896),   south(24, 37, 30.3),    2635.43},
        // ESO VLT, UT4 (ESO website, VLT Whitebook Chapter 2)
        {"VLT4",       "ESO VLT, Paranal, Chile: UT4",            west(70, 24, 8.000),   south(24, 37, 31.0),    2635.43},
        // Gemini South, Cerro Pachon (GPS readings by Patrick Wallace)
        {"GEMINIS",    "Gemini South 8-m telescope",              west(70, 44, 11.5),    south(30, 14, 26.7),    2738.0},
        // Cologne Observatory for Submillimeter Astronomy (KOSMA) (Holger Jakob)
        {"KOSMA3M",    "KOSMA 3m telescope, Gornergrat",          east(7, 47, 3.48),     north(45, 58, 59.772),  3141.0},
        // Magellan 1, 6.5m telescope at Las Campanas, Chile (Skip Schaller)
        {"MAGELLAN1",  "Magellan 1, 6.5m, Las Campanas",          west(70, 41, 31.9),    south(29, 0, 51.7),     2408.0},
        // Magellan 2, 6.5m telescope at Las Campanas, Chile (Skip Schaller)
        {"MAGELLAN2",  "Magellan 2, 6.5m, Las Campanas",          west(70, 41, 33.5),    south(29, 0, 50.3),     2408.0},
        // APEX - Atacama Pathfinder EXperiment, Llano de Chajnantor (APEX web site)
        {"APEX",       "APEX 12m telescope, Llano de Chajnantor", west(67, 45, 33.0),    south(23, 0, 20.8),     5105.0},
        // NANTEN2 Submillimeter Observatory, 4m telescope Atacama desert (NANTEN2 web site)
        {"NANTEN2",    "NANTEN2 4m telescope, Pampa la Bola",     west(67, 42, 08.0),    south(22, 57, 47.0),    4865.0}
    };

    // search for specified observatory
    const Observatory* o = nullptr;
    if (n >= 0) {
        if (n < N_OBS) {
            o = &observatories[n];
        }
    } else if (id != nullptr) {
        size_t nlen = std::strlen(id);
        for (int i = 0; i < N_OBS; i++) {
            const Observatory* p = &observatories[i];
            // if partial `id` is specified, first matching observatory will be returned
            if (std::strncmp(id, p->o_id, nlen) == 0) {
                o = p;
                break;
            }
        }
    }

    // if a match was found, return observatory's data
    if (o != nullptr) {
        obs.o_id = o->o_id;
        obs.o_name = o->o_name;
        obs.o_long = o->o_long;
        obs.o_lat = o->o_lat;
        obs.o_height = o->o_height;
        return true;
    } else {
        obs.o_id = obs.o_name = "?";
        obs.o_long = obs.o_lat = obs.o_height = 0.0;
        return false;
    }
}

}
