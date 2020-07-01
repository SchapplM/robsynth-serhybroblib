% Return the minimum parameter vector for
% palh2m1DE
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
% MPV [21x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 17:39
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = palh2m1DE_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1DE_convert_par2_MPV_fixb: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1DE_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'palh2m1DE_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'palh2m1DE_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t2 = m(5) + m(6);
t1 = m(4) + t2;
t4 = pkin(2) ^ 2;
t9 = t1 * t4;
t3 = pkin(3) ^ 2;
t8 = t2 * t3;
t7 = t3 + t4;
t6 = -mrSges(4,3) - mrSges(5,3);
t5 = pkin(1) ^ 2;
t10 = [(pkin(4) ^ 2 + t5 + t7) * m(6) + (m(3) + m(4) + m(5)) * t5 + t7 * m(5) + t4 * m(4) + Ifges(4,2) + Ifges(5,2) + Ifges(2,3) + Ifges(3,2); (m(3) + t1) * pkin(1) + mrSges(2,1); mrSges(2,2) + mrSges(3,3) - t6; Ifges(3,1) - Ifges(3,2) - t9; Ifges(3,4); t6 * pkin(2) + Ifges(3,5); Ifges(3,6); Ifges(3,3) + t9; t1 * pkin(2) + mrSges(3,1); mrSges(3,2); Ifges(4,1) - Ifges(4,2) - t8; Ifges(4,4); -pkin(3) * mrSges(5,3) + Ifges(4,5); Ifges(4,6); Ifges(4,3) + t8; t2 * pkin(3) + mrSges(4,1); mrSges(4,2); pkin(4) * m(6) + mrSges(5,1); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t10;
