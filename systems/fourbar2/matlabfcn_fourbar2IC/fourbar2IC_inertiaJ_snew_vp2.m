% Calculate joint inertia matrix with Newton Euler for
% fourbar2IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2]';
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
% Datum: 2020-04-24 20:37
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = fourbar2IC_inertiaJ_snew_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(2,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbar2IC_inertiaJ_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'fourbar2IC_inertiaJ_snew_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar2IC_inertiaJ_snew_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'fourbar2IC_inertiaJ_snew_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'fourbar2IC_inertiaJ_snew_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_snew_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:37:51
% EndTime: 2020-04-24 20:37:51
% DurationCPUTime: 0.05s
% Computational Cost: add. (32->12), mult. (28->18), div. (7->3), fcn. (17->4), ass. (0->8)
t44 = cos(qJ(2)) * mrSges(3,1);
t42 = -qJ(3) + qJ(1);
t37 = sin(qJ(2) + t42);
t43 = (-pkin(1) * t37 + pkin(2) * sin(t42)) / t37;
t41 = 0.1e1 / pkin(1);
t39 = sin(qJ(2));
t38 = pkin(2) * mrSges(3,2) * t39;
t1 = [0.2e1 * t38 + Ifges(2,3) + Ifges(3,3) + t39 ^ 2 / t37 ^ 2 * Ifges(4,3) + (m(3) * pkin(2) - 0.2e1 * t44) * pkin(2) + (0.2e1 * (-pkin(2) * t44 + Ifges(3,3) + t38) * t41 + Ifges(3,3) * t41 ^ 2 * t43) * t43;];
Mq = t1;
