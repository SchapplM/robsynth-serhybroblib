% Return the minimum parameter vector for
% fourbar1turnTE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% m [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MPV [25x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 16:23
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = fourbar1turnTE_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnTE_convert_par2_MPV_fixb: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnTE_convert_par2_MPV_fixb: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fourbar1turnTE_convert_par2_MPV_fixb: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'fourbar1turnTE_convert_par2_MPV_fixb: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t1 = pkin(2) ^ 2 * m(4);
t3 = -Ifges(3,2) - t1;
t2 = [pkin(1) ^ 2 * m(5) + Ifges(2,3) - t3; pkin(1) * m(5) + mrSges(2,1); mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3); Ifges(3,1) + t3; Ifges(3,4); -pkin(2) * mrSges(4,3) + Ifges(3,5); Ifges(3,6); Ifges(3,3) + t1; pkin(2) * m(4) + mrSges(3,1); mrSges(3,2); Ifges(4,1) + Ifges(5,2); Ifges(4,4); Ifges(4,5); Ifges(4,2) + Ifges(5,2); Ifges(4,6); Ifges(4,3); mrSges(4,1); mrSges(4,2); Ifges(5,1) - Ifges(5,2); Ifges(5,4); Ifges(5,5); Ifges(5,6); Ifges(5,3); mrSges(5,1); mrSges(5,2);];
MPV = t2;
