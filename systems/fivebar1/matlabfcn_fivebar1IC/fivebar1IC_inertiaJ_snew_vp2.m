% Calculate joint inertia matrix with Newton Euler for
% fivebar1IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AE,BC,CD,ED]';
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
% Datum: 2020-04-27 06:19
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = fivebar1IC_inertiaJ_snew_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fivebar1IC_inertiaJ_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1IC_inertiaJ_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fivebar1IC_inertiaJ_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fivebar1IC_inertiaJ_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'fivebar1IC_inertiaJ_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_snew_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 06:19:11
% EndTime: 2020-04-27 06:19:11
% DurationCPUTime: 0.07s
% Computational Cost: add. (172->37), mult. (234->55), div. (34->7), fcn. (96->11), ass. (0->30)
t92 = sin(qJ(2));
t115 = pkin(2) * t92;
t91 = sin(qJ(4));
t114 = pkin(3) * t91;
t113 = Ifges(3,3) / pkin(5) ^ 2;
t107 = qJ(4) + qJ(3);
t108 = qJ(1) + qJ(2);
t110 = cos(0.2e1 * t107) - cos(0.2e1 * t108);
t81 = 0.1e1 / t110;
t95 = 0.2e1 * qJ(3);
t96 = 0.2e1 * qJ(1);
t112 = (t110 * pkin(4) + (cos(qJ(4) + t95) - cos(qJ(4) + t96 + 0.2e1 * qJ(2))) * pkin(3)) * t81;
t111 = (t110 * pkin(5) + (-cos(qJ(2) + 0.2e1 * qJ(4) + t95) + cos(qJ(2) + t96)) * pkin(2)) * t81;
t109 = Ifges(5,3) / pkin(4) ^ 2;
t106 = mrSges(5,2) * t114;
t105 = cos(qJ(2)) * mrSges(3,1) * pkin(2);
t99 = 0.1e1 / pkin(4);
t104 = t99 * t112;
t97 = 0.1e1 / pkin(5);
t103 = t97 * t111;
t90 = mrSges(3,2) * t115;
t83 = Ifges(3,3) + t90 - t105;
t89 = cos(qJ(4)) * mrSges(5,1) * pkin(3);
t82 = Ifges(5,3) + t89 - t106;
t88 = -sin(t107 - t108);
t85 = 0.1e1 / t88 ^ 2;
t84 = 0.1e1 / t88;
t78 = -Ifges(3,3) * t103 + t83;
t77 = -Ifges(5,3) * t104 + t82;
t1 = [-0.2e1 * t105 + Ifges(2,3) + Ifges(3,3) + 0.2e1 * t90 - (t78 + t83) * t103 + (t85 * t92 ^ 2 * t109 + m(3)) * pkin(2) ^ 2, (t78 * t97 * t114 + (-t109 * t112 + t82 * t99) * t115) * t84; (t77 * t99 * t115 + (-t111 * t113 + t83 * t97) * t114) * t84, -0.2e1 * t106 + Ifges(4,3) + Ifges(5,3) + 0.2e1 * t89 - (t77 + t82) * t104 + (t85 * t91 ^ 2 * t113 + m(5)) * pkin(3) ^ 2;];
Mq = t1;
