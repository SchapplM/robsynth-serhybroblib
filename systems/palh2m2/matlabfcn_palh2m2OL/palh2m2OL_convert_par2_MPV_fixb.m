% Return the minimum parameter vector for
% palh2m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% m [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MPV [38x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 18:09
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = palh2m2OL_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2OL_convert_par2_MPV_fixb: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'palh2m2OL_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'palh2m2OL_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'palh2m2OL_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t101 = m(6) + m(7);
t111 = mrSges(5,3) + mrSges(6,3);
t99 = m(5) + t101;
t97 = m(4) + t99;
t94 = pkin(4) ^ 2 * t97;
t110 = -Ifges(3,2) - t94;
t95 = pkin(2) ^ 2 * t99;
t109 = -Ifges(4,2) - t95;
t98 = pkin(5) ^ 2 * t101;
t108 = -Ifges(5,2) - t98;
t100 = pkin(3) ^ 2 * m(7);
t107 = -Ifges(6,2) - t100;
t106 = mrSges(4,3) + t111;
t96 = m(3) + t97;
t1 = [pkin(1) ^ 2 * t96 + Ifges(2,3) - t107 - t108 - t109 - t110; pkin(1) * t96 + mrSges(2,1); mrSges(2,2) - mrSges(3,3) - t106; Ifges(3,1) + t110; Ifges(3,4); -pkin(4) * t106 + Ifges(3,5); Ifges(3,6); Ifges(3,3) + t94; pkin(4) * t97 + mrSges(3,1); mrSges(3,2); Ifges(4,1) + t109; Ifges(4,4); -pkin(2) * t111 + Ifges(4,5); Ifges(4,6); Ifges(4,3) + t95; pkin(2) * t99 + mrSges(4,1); mrSges(4,2); Ifges(5,1) + t108; Ifges(5,4); -pkin(5) * mrSges(6,3) + Ifges(5,5); Ifges(5,6); Ifges(5,3) + t98; pkin(5) * t101 + mrSges(5,1); mrSges(5,2); Ifges(6,1) + Ifges(7,2) + t107; -pkin(3) * mrSges(7,3) + Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(7,2) + Ifges(6,3) + t100; pkin(3) * m(7) + mrSges(6,1); mrSges(6,2) + mrSges(7,3); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV = t1;
