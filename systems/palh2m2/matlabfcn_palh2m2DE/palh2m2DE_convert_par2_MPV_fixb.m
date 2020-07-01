% Return the minimum parameter vector for
% palh2m2DE
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
% MPV [22x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 17:56
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = palh2m2DE_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2DE_convert_par2_MPV_fixb: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'palh2m2DE_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'palh2m2DE_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'palh2m2DE_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t3 = m(6) + m(7);
t8 = m(5) + t3;
t2 = m(4) + t8;
t6 = pkin(4) ^ 2;
t11 = t2 * t6;
t5 = pkin(5) ^ 2;
t10 = t3 * t5;
t4 = pkin(2) ^ 2;
t9 = t4 + t5;
t7 = -mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t1 = m(3) + t2;
t12 = [t1 * pkin(1) ^ 2 + (pkin(3) ^ 2 + t6 + t9) * m(7) + (m(4) + m(5) + m(6)) * t6 + t9 * m(6) + t4 * m(5) + Ifges(5,2) + Ifges(6,2) + Ifges(2,3) + Ifges(3,2) + Ifges(4,2); t1 * pkin(1) + mrSges(2,1); mrSges(2,2) - mrSges(3,3) + t7; Ifges(3,1) - Ifges(3,2) - t11; Ifges(3,4); t7 * pkin(4) + Ifges(3,5); Ifges(3,6); Ifges(3,3) + t11; t2 * pkin(4) + mrSges(3,1); mrSges(3,2); t8 * pkin(2) + mrSges(4,1); Ifges(5,1) - Ifges(5,2) - t10; Ifges(5,4); -pkin(5) * mrSges(6,3) + Ifges(5,5); Ifges(5,6); Ifges(5,3) + t10; t3 * pkin(5) + mrSges(5,1); mrSges(5,2); pkin(3) * m(7) + mrSges(6,1); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV = t12;
