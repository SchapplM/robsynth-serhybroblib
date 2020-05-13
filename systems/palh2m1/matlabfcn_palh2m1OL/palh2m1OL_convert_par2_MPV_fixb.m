% Return the minimum parameter vector for
% palh2m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MPV [31x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 00:53
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = palh2m1OL_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1OL_convert_par2_MPV_fixb: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1OL_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'palh2m1OL_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'palh2m1OL_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t96 = m(5) + m(6);
t105 = mrSges(4,3) + mrSges(5,3);
t94 = m(4) + t96;
t91 = pkin(2) ^ 2 * t94;
t104 = (-Ifges(3,2) - t91);
t93 = pkin(3) ^ 2 * t96;
t103 = (-Ifges(4,2) - t93);
t102 = 2 * pkin(6) * mrSges(6,3) + Ifges(6,2);
t101 = m(6) * pkin(6) + mrSges(6,3);
t98 = (pkin(4) ^ 2);
t97 = pkin(6) ^ 2;
t92 = (m(3) + t94);
t1 = [pkin(1) ^ 2 * t92 + t98 * m(6) + Ifges(5,2) + Ifges(2,3) - t103 - t104; pkin(1) * t92 + mrSges(2,1); mrSges(2,2) + mrSges(3,3) + t105; Ifges(3,1) + t104; Ifges(3,4); -pkin(2) * t105 + Ifges(3,5); Ifges(3,6); Ifges(3,3) + t91; pkin(2) * t94 + mrSges(3,1); mrSges(3,2); Ifges(4,1) + t103; Ifges(4,4); -mrSges(5,3) * pkin(3) + Ifges(4,5); Ifges(4,6); Ifges(4,3) + t93; pkin(3) * t96 + mrSges(4,1); mrSges(4,2); Ifges(5,1) - Ifges(5,2) + (t97 - t98) * m(6) + t102; t101 * pkin(4) + Ifges(5,4); Ifges(5,5); Ifges(5,6); Ifges(5,3) + (t97 + t98) * m(6) + t102; pkin(4) * m(6) + mrSges(5,1); mrSges(5,2) - t101; Ifges(6,1) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t1;
