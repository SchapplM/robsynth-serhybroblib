% Calculate joint inertia matrix with Newton Euler for
% fourbar1turnIC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
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
% Mq [2x2]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 11:33
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = fourbar1turnIC_inertiaJ_snew_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnIC_inertiaJ_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnIC_inertiaJ_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnIC_inertiaJ_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fourbar1turnIC_inertiaJ_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'fourbar1turnIC_inertiaJ_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_snew_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 11:33:36
% EndTime: 2020-05-07 11:33:37
% DurationCPUTime: 0.14s
% Computational Cost: add. (119->57), mult. (153->83), div. (16->5), fcn. (142->8), ass. (0->24)
t100 = sin(qJ(3));
t103 = cos(qJ(3));
t115 = -t103 * mrSges(4,1) + t100 * mrSges(4,2);
t114 = pkin(2) * t100;
t113 = pkin(2) * t103;
t104 = cos(qJ(2));
t112 = pkin(2) * t104;
t101 = sin(qJ(2));
t92 = t100 * t101 - t103 * t104;
t93 = -t100 * t104 - t103 * t101;
t111 = Ifges(4,5) * t93 + Ifges(4,6) * t92;
t108 = -qJ(4) + qJ(2);
t97 = sin(qJ(3) + t108);
t96 = 0.1e1 / t97;
t107 = (-pkin(3) * t97 + pkin(2) * sin(t108)) / pkin(3) * t96;
t106 = pkin(2) * t115 + Ifges(4,3);
t102 = cos(qJ(4));
t99 = sin(qJ(4));
t95 = mrSges(4,3) * t113 + Ifges(4,5);
t94 = -mrSges(4,3) * t114 + Ifges(4,6);
t88 = 0.1e1 / pkin(4) * t96 * (t99 * Ifges(5,5) + t102 * Ifges(5,6)) * t114;
t87 = -mrSges(4,2) * t112 + Ifges(4,1) * t93 + Ifges(4,4) * t92;
t86 = mrSges(4,1) * t112 + Ifges(4,4) * t93 + Ifges(4,2) * t92;
t1 = [Ifges(2,3) + t104 * (Ifges(3,2) * t104 - t100 * t87 - t103 * t86 - pkin(2) * (-m(4) * t112 - t92 * mrSges(4,1) + t93 * mrSges(4,2))) + Ifges(5,1) * t99 ^ 2 + (0.2e1 * Ifges(5,4) * t99 + Ifges(5,2) * t102) * t102 + (Ifges(3,1) * t101 + 0.2e1 * Ifges(3,4) * t104 + t100 * t86 - t103 * t87) * t101 + (m(5) * pkin(1) + 0.2e1 * t102 * mrSges(5,1) - 0.2e1 * t99 * mrSges(5,2)) * pkin(1), t101 * (t100 * t94 - t103 * t95 + Ifges(3,5)) + t104 * (-t100 * t95 - t103 * t94 + Ifges(3,6)) + (t101 * (-t103 * Ifges(4,5) + t100 * Ifges(4,6)) + t104 * (-t100 * Ifges(4,5) - t103 * Ifges(4,6))) * t107 + t88; Ifges(3,5) * t101 + Ifges(3,6) * t104 + t111 * t107 + t88 + pkin(2) * (-t100 * t92 + t103 * t93) * mrSges(4,3) + t111, Ifges(3,3) + t106 + (Ifges(4,3) * t107 + Ifges(4,3) + t106) * t107 + (-t103 * (-m(4) * t113 + mrSges(4,1)) + t115 * t107 + (mrSges(4,2) + (m(4) + 0.1e1 / pkin(4) ^ 2 / t97 ^ 2 * Ifges(5,3)) * t114) * t100) * pkin(2);];
Mq = t1;
