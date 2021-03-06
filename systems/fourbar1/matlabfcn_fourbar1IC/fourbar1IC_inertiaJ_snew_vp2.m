% Calculate joint inertia matrix with Newton Euler for
% fourbar1IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4]';
% m [4x1]
%   mass of all robot links (including the base)
% mrSges [4x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [1x1]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 20:15
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = fourbar1IC_inertiaJ_snew_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbar1IC_inertiaJ_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1IC_inertiaJ_snew_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar1IC_inertiaJ_snew_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'fourbar1IC_inertiaJ_snew_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'fourbar1IC_inertiaJ_snew_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_snew_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:15:51
% EndTime: 2020-04-24 20:15:51
% DurationCPUTime: 0.05s
% Computational Cost: add. (32->12), mult. (31->19), div. (8->4), fcn. (17->4), ass. (0->8)
t44 = -qJ(3) + qJ(1);
t37 = sin(qJ(2) + t44);
t45 = (-pkin(3) * t37 + pkin(2) * sin(t44)) / t37;
t43 = mrSges(3,1) * cos(qJ(2)) * pkin(2);
t41 = 0.1e1 / pkin(3);
t39 = sin(qJ(2));
t38 = pkin(2) * mrSges(3,2) * t39;
t1 = [-0.2e1 * t43 + Ifges(2,3) + Ifges(3,3) + 0.2e1 * t38 + (m(3) + t39 ^ 2 / pkin(4) ^ 2 / t37 ^ 2 * Ifges(4,3)) * pkin(2) ^ 2 + (0.2e1 * (Ifges(3,3) + t38 - t43) * t41 + Ifges(3,3) * t41 ^ 2 * t45) * t45;];
Mq = t1;
